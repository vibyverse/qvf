#!/usr/bin/env python3
"""Typed C++ payload/kind rendering used by gen_kinds.py (Apache-2.0).

Schema assertions remain in the shared generated checker. This module maps
known object fields to C++ scalar/array/record types and preserves open fields
in `extra`. It never substitutes a whole canonical payload with a JSON wrapper.
"""
from __future__ import annotations
import json,re

def ident(name):
    out=re.sub(r'[^A-Za-z0-9_]', '_', name)
    if not out or out[0].isdigit(): out='field_'+out
    if out in {'class','struct','template','namespace','operator','default','delete','new','private','public','protected','virtual','typename','union','using','const','volatile','static','inline','explicit','friend','return','while','for','if','else','switch','case','break','continue','do','try','catch','throw','auto','register','extern','mutable','thread_local','alignas','alignof','noexcept','sizeof','typeid','this','true','false','nullptr','enum','bool','char','double','float','int','long','short','signed','unsigned','void','wchar_t','char16_t','char32_t','and','or','xor','not','bitand','bitor','compl','and_eq','or_eq','xor_eq','not_eq','to_json','from_json','to_section','from_section','extra'}: out+='_'
    return out

def pascal(name): return ''.join(x[:1].upper()+x[1:] for x in re.split(r'[^A-Za-z0-9]+',name) if x)
def quote(s): return json.dumps(s,ensure_ascii=False)

class Types:
    def __init__(self,schema): self.defs=schema['$defs'];self.done={};self.lines=[];self.fallbacks=[]
    def merged(self,node):
        if not isinstance(node,dict): return {}
        result={}; props={};required=[]
        parts=[]
        if '$ref' in node: parts.append(self.defs[node['$ref'].rsplit('/',1)[-1]])
        parts+=node.get('allOf',[])
        for p in parts+[node]:
            p=self.merged(p) if p is not node else p
            for k,v in p.items():
                if k not in ('$ref','allOf','properties','required'): result[k]=v
            for k,v in p.get('properties',{}).items():
                if k in props: props[k]={'allOf':[props[k],v]}
                else: props[k]=v
            required+=p.get('required',[])
        if props: result['properties']=props
        if required: result['required']=list(dict.fromkeys(required))
        return result
    def union(self,types):
        types=list(dict.fromkeys(types))
        if len(types)==1: return types[0]
        if 'std::monostate' in types and len(types)==2:
            return 'std::optional<'+next(t for t in types if t!='std::monostate')+'>'
        return 'std::variant<'+', '.join(types)+'>'
    def type(self,node,name):
        if not isinstance(node,dict) or not node: self.fallbacks.append(name);return 'Json'
        if '$ref' in node and all(k in ('$ref','description','$comment','deprecated') for k in node):
            ref=node['$ref'].rsplit('/',1)[-1]
            return self.named(self.defs[ref],ref)
        node=self.merged(node)
        typ=node.get('type')
        if isinstance(typ,list): return self.union([self.type(dict(node,type=t),name+pascal(t)) for t in typ])
        if not typ and 'const' in node: typ=self.literal_type(node['const'])
        if not typ and node.get('enum'): typ=self.literal_type(node['enum'][0])
        if not typ and node.get('properties'): typ='object'
        if not typ and ('anyOf' in node or 'oneOf' in node):
            choices=node.get('anyOf',node.get('oneOf'))
            return self.union([self.type(n,name+'Alternative'+str(i+1)) for i,n in enumerate(choices)])
        basic={'string':'std::string','integer':'int64_t','number':'double','boolean':'bool','null':'std::monostate'}
        if typ in basic:return basic[typ]
        if typ=='array':
            if node.get('prefixItems'):
                types=[self.type(x,name+'Item'+str(i+1)) for i,x in enumerate(node['prefixItems'])]
                if len(set(types))==1 and node.get('maxItems')==len(types):return f'std::array<{types[0]}, {len(types)}>'
                if node.get('items') is False or node.get('maxItems')==len(types):return 'std::tuple<'+', '.join(types)+'>'
                # Heterogeneous prefixes followed by a homogeneous tail are
                # represented by a vector of typed alternatives.
                item=self.union(types+[self.type(node.get('items',{}),name+'Item')])
            else:item=self.type(node.get('items',{}),name+'Item')
            if node.get('minItems')==node.get('maxItems') and 'maxItems' in node:
                return f'std::array<{item}, {node["maxItems"]}>'
            return f'std::vector<{item}>'
        if typ=='object' or 'patternProperties' in node:
            props=node.get('properties',{})
            if not props:
                subs=list(node.get('patternProperties',{}).values())
                ap=node.get('additionalProperties')
                if isinstance(ap,dict):subs.append(ap)
                item=self.union([self.type(n,name+'Entry'+str(i)) for i,n in enumerate(subs)]) if subs else 'Json'
                return f'std::map<std::string, {item}>'
            return self.record(node,name)
        self.fallbacks.append(name);return 'Json'
    @staticmethod
    def literal_type(x):
        if x is None:return 'null'
        if isinstance(x,bool):return 'boolean'
        if isinstance(x,str):return 'string'
        if isinstance(x,int):return 'integer'
        if isinstance(x,float):return 'number'
        raise ValueError('non-scalar schema enum/const needs a typed representation')
    def named(self,node,name):
        if name in self.done:return self.done[name]
        t=self.type(node,name)
        self.done[name]=t
        return t
    def record(self,node,name):
        if name in self.done:return self.done[name]
        self.done[name]=name
        props=node.get('properties',{}); req=set(node.get('required',[]))
        fields=[(k,self.type(v,name+pascal(k)),k in req) for k,v in props.items()]
        L=['struct '+name+' {']
        for k,t,required in fields:
            default='{}'
            merged=self.merged(props[k])
            if required and 'const' in merged and isinstance(merged['const'],(str,bool,int,float)):
                default=' = '+quote(merged['const'])
            L.append(f'    {t if required else "std::optional<"+t+">"} {ident(k)}{default};')
        L+=['    std::map<std::string, Json> extra;', '    Json to_json() const {','        Json out = codec::Value<std::map<std::string, Json>>::encode(extra);']
        for k,t,required in fields:
            L.append(f'        if (out.has({quote(k)})) throw std::invalid_argument("use the typed field, not extra: {k}");')
            L.append(f'        codec::put(out, {quote(k)}, {ident(k)});' if required else f'        if ({ident(k)}) codec::put(out, {quote(k)}, *{ident(k)});')
        L+=['        return out;','    }',f'    static {name} from_json(const Json& j) {{','        codec::require(j.is_object(), "expected typed payload object");',f'        {name} out;']
        for k,t,required in fields:
            L.append(f'        out.{ident(k)} = codec::get<{t}>(j, {quote(k)});' if required else f'        if (j.has({quote(k)})) out.{ident(k)}.emplace(codec::get<{t}>(j, {quote(k)}));')
        checks=' && '.join(f'kv.first != {quote(k)}' for k in props) or 'true'
        L+=[f'        for (const auto& kv : j.items()) if ({checks}) out.extra[kv.first] = kv.second;', '        return out;','    }','};','']
        self.lines += L
        return name

NAMES={'wavefunction.gto':'WavefunctionGTO','wavefunction.bloch':'WavefunctionBloch','job.spec':'JobSpec','run.record':'RunRecord','dos.total':'DosTotal','dos.projected':'DosProjected','dos.coop':'DosCoop','dos.cohp':'DosCohp','topology.qtaim':'TopologyQTAIM','spectra.ir':'SpectraIR','spectra.raman':'SpectraRaman','spectra.uvvis':'SpectraUVVis','spectra.ecd':'SpectraECD','spectra.vcd':'SpectraVCD','spectra.nmr':'SpectraNMR','spectra.epr':'SpectraEPR','basis.ao':'BasisAO','volume.elf':'VolumeELF','volume.rdg':'VolumeRDG','scf_history':'ScfHistory'}
def render(schema,model):
    types=Types(schema); sections=[]
    common=['id','kind','label','component','critical','structure_ref','run_ref','group','members']
    for kind in model.kinds:
        name=NAMES.get(kind['kind'],pascal(kind['kind']))
        node=types.merged(schema['$defs'][kind['def']]); props={k:v for k,v in node.get('properties',{}).items() if k not in ('kind','members')}
        peers=types.record({'type':'object','properties':props,'required':node.get('required',[])},name+'Fields') if props else None
        # Each kind has distinct named members and payload types, even when two
        # kinds use the same payload schema (volumes, spectra, DOS).
        roles=[]
        for role in kind['roles']:
            if role['format']=='json': t=types.named(schema['$defs'][role['payload_def']],role['payload_def']) if role['payload_def'] else 'Json'
            elif role['format']=='binary': t='TensorData'
            else:t='OpaqueData'
            field=ident(role['role']) if not role['pattern'] else ('attachments' if 'attachment' in role['role'] else 'custom_members')
            if role['pattern']: t='std::map<std::string, '+t+'>'
            roles.append((role,field,t))
        L=[f'struct {name} : SectionFields {{', f'    static constexpr std::string_view kind = {quote(kind["kind"])};',f'    {name}() {{ id = {quote(kind["kind"].replace(".","_"))}; }}']
        if peers:L.append(f'    {peers} fields;')
        for r,f,t in roles:L.append(f'    {t if r["required"] or r["pattern"] else "std::optional<"+t+">"} {f}{{}};')
        L+=['    std::string to_section(Writer& writer) const {','        Section s = descriptor(std::string(kind));']
        if peers:L.append('        const Json peer_fields = fields.to_json();\n        for (const auto& kv : peer_fields.items()) s.extra.set(kv.first, kv.second);')
        L+=['        std::vector<JsonPayload> json; std::vector<BinaryPayload> binary; std::vector<OpaquePayload> opaque;']
        for r,f,t in roles:
            container={'json':'json','binary':'binary','opaque':'opaque'}[r['format']]
            def push(role,value):
                if r['format']=='json':return f'json.push_back({{{role}, codec::Value<{t if not r["pattern"] else t[len("std::map<std::string, "):-1]}>::encode({value})}});'
                if r['format']=='binary':return f'binary.push_back({{{role}, {value}}});'
                return f'opaque.push_back({{{role}, {value}.bytes, {value}.media_type}});'
            if r['pattern']:L.append(f'        for (const auto& member : {f}) '+push('member.first','member.second'))
            elif r['required']:L.append('        '+push(quote(r['role']),f))
            else:L.append(f'        if ({f}) '+push(quote(r['role']),'(*'+f+')'))
        L+=['        return writer.add_payload(std::move(s), json, binary, opaque);','    }',f'    static {name} from_section(const Reader& reader, const Section& section) {{','        if (!section.kind || *section.kind != kind) throw std::invalid_argument("wrong section kind");',f'        {name} out; out.load(section);']
        if peers:
            L.append(f'        out.fields = {peers}::from_json(section.extra);')
            L.append('        out.fields.extra.clear();')
            L.append('        out.extra = Json::object();')
            checks=' && '.join(f'kv.first != {quote(k)}' for k in props) or 'true'
            L.append(f'        for (const auto& kv : section.extra.items()) if ({checks}) out.extra.set(kv.first, kv.second);')
        for r,f,t in roles:
            def read(role):
                if r['format']=='json':return f'codec::Value<{t if not r["pattern"] else t[len("std::map<std::string, "):-1]}>::decode(reader.json(section, {role}))'
                if r['format']=='binary':return f'TensorData(reader.tensor(section, {role}))'
                return f'OpaqueData{{reader.bytes(section, {role}), section.member({role})->media_type}}'
            if r['pattern']:
                # Pattern identities are generated; the shared schema checks
                # exact regex validity when writing the reconstituted section.
                known=[rr['role'] for rr,_,_ in roles if not rr['pattern']]
                checks=' && '.join(f'member.first != {quote(x)}' for x in known) or 'true'
                L.append(f'        if (section.members) for (const auto& member : *section.members) if ({checks}) out.{f}.emplace(member.first, '+read('member.first')+');')
            elif r['required']:L.append(f'        out.{f} = '+read(quote(r['role']))+';')
            else:L.append(f'        if (section.member({quote(r["role"])})) out.{f}.emplace('+read(quote(r['role']))+');')
        L+=['        return out;','    }','};','']
        sections+=L
    head=['// GENERATED by scripts/gen_kinds.py from spec/qvf_manifest.schema.json; do not edit','// License: Apache-2.0.','#ifndef QVF_KINDS_ALL_HPP','#define QVF_KINDS_ALL_HPP','#include "qvf/kinds/codec.hpp"','#include "qvf/reader.hpp"','#include "qvf/writer.hpp"','namespace qvf { namespace kinds {','struct OpaqueData { std::vector<uint8_t> bytes; std::optional<std::string> media_type; };','struct SectionFields {','    std::string id;','    std::optional<std::string> label, component, structure_ref, run_ref, group;','    std::optional<bool> critical;','    Json extra = Json::object();','    Section descriptor(const std::string& kind) const {','        Section s; s.id=id; s.kind=kind; s.label=label; s.component=component;','        s.structure_ref=structure_ref; s.run_ref=run_ref; s.group=group; s.critical=critical; s.extra=extra; return s;','    }','    void load(const Section& s) {','        if (!s.id) throw std::invalid_argument("section id absent");','        id=*s.id; label=s.label; component=s.component; structure_ref=s.structure_ref; run_ref=s.run_ref; group=s.group; critical=s.critical; extra=s.extra;','    }','};','']
    copy=['inline std::string copy_section(const Reader& reader, const Section& section, Writer& writer) {', '    if (!section.kind || !section.id) throw std::invalid_argument("section id/kind absent");']
    for kind in model.kinds:
        name=NAMES.get(kind['kind'],pascal(kind['kind']))
        copy.append(f'    if (*section.kind == {quote(kind["kind"])}) return writer.add({name}::from_section(reader, section));')
    copy += ['    throw std::invalid_argument("copy_section requires a canonical kind");', '}']
    text='\n'.join(head+types.lines+sections+copy+['} } // qvf::kinds','#endif',''])
    return text,types.fallbacks


# --- Rust: rust/src/kinds.rs ---------------------------------------------------
#
# The same mapping as the C++ records, as serde models: known fields get Rust
# types, open fields land in a flattened `extra` map, and every canonical kind
# gets a struct with its section fields and typed members that implements
# `crate::typed::TypedKind`.

RUST_KEYWORDS = {'as', 'async', 'await', 'break', 'const', 'continue', 'crate', 'dyn', 'else', 'enum', 'extern',
                 'false', 'fn', 'for', 'if', 'impl', 'in', 'let', 'loop', 'match', 'mod', 'move', 'mut', 'pub',
                 'ref', 'return', 'self', 'Self', 'static', 'struct', 'super', 'trait', 'true', 'type', 'unsafe',
                 'use', 'where', 'while', 'abstract', 'become', 'box', 'do', 'final', 'macro', 'override', 'priv',
                 'typeof', 'unsized', 'virtual', 'yield', 'try', 'gen', 'extra'}


def rust_ident(name):
    out = re.sub(r'[^A-Za-z0-9_]', '_', name)
    if not out or out[0].isdigit():
        out = 'field_' + out
    if out in RUST_KEYWORDS:
        out += '_'
    return out


def rust_quote(s):
    return json.dumps(s, ensure_ascii=False)


class RustTypes(Types):
    """Rust type names for schema nodes, collecting struct and enum items."""

    def __init__(self, schema):
        super().__init__(schema)
        self.unions = {}
        self.used = set()
        self.stack = []

    def unique(self, name):
        base, n = name, 2
        while name in self.used:
            name = f'{base}{n}'
            n += 1
        self.used.add(name)
        return name

    @staticmethod
    def variant(t):
        simple = {'String': 'Str', 'i64': 'Int', 'f64': 'Num', 'bool': 'Bool', '()': 'Null', 'Value': 'Json'}
        if t in simple:
            return simple[t]
        for prefix, v in (('Vec<', 'List'), ('BTreeMap<', 'Map'), ('[', 'Array'), ('(', 'Tuple'), ('Option<', 'Maybe')):
            if t.startswith(prefix):
                return v
        return t

    def union(self, types, name='Choice'):
        types = list(dict.fromkeys(types))
        if len(types) == 1:
            return types[0]
        if '()' in types and len(types) == 2:
            return 'Option<' + next(t for t in types if t != '()') + '>'
        # An untagged enum tries variants in order: integers before numbers,
        # so that an integer stays an integer through a round trip.
        if 'i64' in types and 'f64' in types:
            types.remove('i64')
            types.insert(types.index('f64'), 'i64')
        key = tuple(types)
        if key in self.unions:
            return self.unions[key]
        enum = self.unique(name)
        self.unions[key] = enum
        seen = set()
        L = ['#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]', '#[serde(untagged)]', f'pub enum {enum} {{']
        for t in types:
            v = self.variant(t)
            base, n = v, 2
            while v in seen:
                v = f'{base}{n}'
                n += 1
            seen.add(v)
            L.append(f'    {v},' if t == '()' else f'    {v}({self.boxed(t)}),')
        L += ['}', '']
        self.lines += L
        return enum

    def boxed(self, t):
        # A record that contains itself (directly or through an enum) needs
        # indirection; Vec and BTreeMap already provide it.
        return f'Box<{t}>' if t in self.stack else t

    def type(self, node, name):
        if not isinstance(node, dict) or not node:
            self.fallbacks.append(name)
            return 'Value'
        if '$ref' in node and all(k in ('$ref', 'description', '$comment', 'deprecated') for k in node):
            ref = node['$ref'].rsplit('/', 1)[-1]
            return self.named(self.defs[ref], ref)
        node = self.merged(node)
        typ = node.get('type')
        if isinstance(typ, list):
            return self.union([self.type(dict(node, type=t), name + pascal(t)) for t in typ], name)
        if not typ and 'const' in node:
            typ = self.literal_type(node['const'])
        if not typ and node.get('enum'):
            typ = self.literal_type(node['enum'][0])
        if not typ and node.get('properties'):
            typ = 'object'
        if not typ and ('anyOf' in node or 'oneOf' in node):
            choices = node.get('anyOf', node.get('oneOf'))
            return self.union([self.type(n, name + 'Alternative' + str(i + 1)) for i, n in enumerate(choices)], name)
        basic = {'string': 'String', 'integer': 'i64', 'number': 'f64', 'boolean': 'bool', 'null': '()'}
        if typ in basic:
            return basic[typ]
        if typ == 'array':
            if node.get('prefixItems'):
                types = [self.type(x, name + 'Item' + str(i + 1)) for i, x in enumerate(node['prefixItems'])]
                if len(set(types)) == 1 and node.get('maxItems') == len(types) and len(types) <= 32:
                    return f'[{types[0]}; {len(types)}]'
                if (node.get('items') is False or node.get('maxItems') == len(types)) and len(types) <= 16:
                    return '(' + ', '.join(types) + ',)' if len(types) == 1 else '(' + ', '.join(types) + ')'
                item = self.union(types + [self.type(node.get('items', {}), name + 'Item')], name + 'Item')
            else:
                item = self.type(node.get('items', {}), name + 'Item')
            if 'maxItems' in node and node.get('minItems') == node['maxItems'] and 0 < node['maxItems'] <= 32:
                return f'[{item}; {node["maxItems"]}]'
            return f'Vec<{item}>'
        if typ == 'object' or 'patternProperties' in node:
            props = node.get('properties', {})
            if not props:
                subs = list(node.get('patternProperties', {}).values())
                ap = node.get('additionalProperties')
                if isinstance(ap, dict):
                    subs.append(ap)
                item = self.union([self.type(n, name + 'Entry' + str(i)) for i, n in enumerate(subs)], name + 'Entry') if subs else 'Value'
                return f'BTreeMap<String, {item}>'
            return self.record(node, name)
        self.fallbacks.append(name)
        return 'Value'

    def record(self, node, name):
        if name in self.done:
            return self.done[name]
        self.used.add(name)
        self.done[name] = name
        self.stack.append(name)
        props = node.get('properties', {})
        req = set(node.get('required', []))
        fields = [(k, self.type(v, name + pascal(k)), k in req) for k, v in props.items()]
        L = ['#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]', f'pub struct {name} {{']
        idents = set()
        for k, t, required in fields:
            ident = rust_ident(k)
            while ident in idents:
                ident += '_'
            idents.add(ident)
            attrs = [] if ident == k else [f'rename = {rust_quote(k)}']
            t = self.boxed(t)
            if not required:
                attrs += ['default', 'skip_serializing_if = "Option::is_none"']
                if t.startswith('Option<'):
                    # Absent and null are different: absent is the outer None.
                    attrs.append('with = "crate::typed::double_option"')
                t = f'Option<{t}>'
            if attrs:
                L.append(f'    #[serde({", ".join(attrs)})]')
            L.append(f'    pub {ident}: {t},')
        L += ['    /// Fields the schema leaves open, kept for the round trip.',
              '    #[serde(flatten)]', '    pub extra: BTreeMap<String, Value>,', '}', '']
        self.stack.pop()
        self.lines += L
        return name


def render_rust(schema, model):
    types = RustTypes(schema)
    for kind in model.kinds:
        types.used.add(NAMES.get(kind['kind'], pascal(kind['kind'])))
    items = []
    names = []
    for kind in model.kinds:
        name = NAMES.get(kind['kind'], pascal(kind['kind']))
        node = types.merged(schema['$defs'][kind['def']])
        # The common section fields of `$defs/Section` (id, label, refs) come
        # first, then the kind's own; `kind` and `members` are the struct.
        base = schema['$defs']['Section']
        props = {k: v for k, v in list(base['properties'].items()) + list(node.get('properties', {}).items())
                 if k not in ('kind', 'members')}
        required = [k for k in base.get('required', []) + node.get('required', []) if k in props]
        fields = types.record({'type': 'object', 'properties': props, 'required': required}, name + 'Fields')
        roles = []
        known = [r['role'] for r in kind['roles'] if not r['pattern']]
        patterns = [r for r in kind['roles'] if r['pattern']]
        if len(patterns) > 1:
            raise ValueError(f'{kind["kind"]}: more than one pattern role has no typed representation')
        for role in kind['roles']:
            if role['format'] == 'json':
                t = types.named(schema['$defs'][role['payload_def']], role['payload_def']) if role['payload_def'] else 'Value'
            elif role['format'] == 'binary':
                t = 'Tensor'
            else:
                t = 'OpaqueData'
            field = rust_ident(role['role']) if not role['pattern'] else ('attachments' if 'attachment' in role['role'] else 'custom_members')
            roles.append((role, field, t))
        doc = f'/// A `{kind["kind"]}` section with typed members (`{kind["def"]}`).'
        L = [doc, '#[derive(Debug, Clone, PartialEq)]', f'pub struct {name} {{',
             '    /// The section fields other than `kind` and `members`.', f'    pub fields: {fields},']
        for r, f, t in roles:
            if r['pattern']:
                L.append(f'    /// Members whose role matches `{r["role"]}`, by role.')
                L.append(f'    pub {f}: BTreeMap<String, {t}>,')
            else:
                L.append(f'    pub {f}: {t if r["required"] else "Option<" + t + ">"},')
        L += ['}', '', f'impl TypedKind for {name} {{', f'    const KIND: &\'static str = {rust_quote(kind["kind"])};', '',
              '    fn from_section(reader: &Reader, section: Section<\'_>) -> Result<Self, QvfError> {',
              '        check_kind(section, Self::KIND)?;']
        reader_fn = {'json': 'json_member', 'binary': 'tensor_member', 'opaque': 'opaque_member'}
        for r, f, t in roles:
            fn = reader_fn[r['format']]
            if r['pattern']:
                quoted = ', '.join(rust_quote(x) for x in known)
                L.append(f'        let {f} = pattern_members(section, &[{quoted}], |role| {fn}(reader, section, role))?;')
            elif r['required']:
                L.append(f'        let {f} = {fn}(reader, section, {rust_quote(r["role"])})?;')
            else:
                L.append(f'        let {f} = optional(section, {rust_quote(r["role"])}, |role| {fn}(reader, section, role))?;')
        L.append('        Ok(Self { fields: section_fields(section)?, ' + ''.join(f'{f}, ' for _, f, _ in roles) + '})')
        L += ['    }', '', '    fn to_writer<W: Write>(&self, writer: &mut Writer<W>) -> Result<String, WriteError> {',
              '        let id = add_section(writer, Self::KIND, &self.fields)?;']
        writer_fn = {'json': 'put_json', 'binary': 'put_tensor', 'opaque': 'put_opaque'}
        for r, f, t in roles:
            fn = writer_fn[r['format']]
            if r['pattern']:
                L.append(f'        for (role, value) in &self.{f} {{')
                L.append(f'            {fn}(writer, &id, role, value)?;')
                L.append('        }')
            elif r['required']:
                L.append(f'        {fn}(writer, &id, {rust_quote(r["role"])}, &self.{f})?;')
            else:
                L.append(f'        if let Some(value) = &self.{f} {{')
                L.append(f'            {fn}(writer, &id, {rust_quote(r["role"])}, value)?;')
                L.append('        }')
        L += ['        Ok(id)', '    }', '', '    fn id(&self) -> &str {', '        &self.fields.id', '    }', '}', '']
        items += L
        names.append((kind['kind'], name))
    any_kind = ['/// A section of any canonical kind, decoded with its typed model.',
                '#[derive(Debug, Clone, PartialEq)]', '#[allow(clippy::large_enum_variant)]', 'pub enum AnyKind {']
    any_kind += [f'    {n}({n}),' for _, n in names]
    any_kind += ['}', '', 'impl AnyKind {',
                 '    /// Decode a section of a canonical kind. Returns `Ok(None)` for a vendor kind.',
                 '    pub fn from_section(reader: &Reader, section: Section<\'_>) -> Result<Option<Self>, QvfError> {',
                 '        Ok(Some(match section.kind {']
    any_kind += [f'            {n}::KIND => AnyKind::{n}({n}::from_section(reader, section)?),' for _, n in names]
    any_kind += ['            _ => return Ok(None),', '        }))', '    }', '',
                 '    /// The canonical kind.', '    pub fn kind(&self) -> &\'static str {', '        match self {']
    any_kind += [f'            AnyKind::{n}(_) => {n}::KIND,' for _, n in names]
    any_kind += ['        }', '    }', '', '    /// The section id.', '    pub fn id(&self) -> &str {', '        match self {']
    any_kind += [f'            AnyKind::{n}(v) => v.id(),' for _, n in names]
    any_kind += ['        }', '    }', '',
                 '    /// Write the section and its members; returns its id.',
                 '    pub fn to_writer<W: Write>(&self, writer: &mut Writer<W>) -> Result<String, WriteError> {',
                 '        match self {']
    any_kind += [f'            AnyKind::{n}(v) => v.to_writer(writer),' for _, n in names]
    any_kind += ['        }', '    }', '}', '']
    head = ['// GENERATED by scripts/gen_kinds.py from spec/qvf_manifest.schema.json; do not edit',
            '// License: Apache-2.0.',
            '//! Typed serde models of the schema `$defs` payloads and of every canonical',
            '//! section kind, generated by `scripts/gen_kinds.py` (the same mapping as',
            '//! `cpp/include/qvf/kinds/all.hpp`). Fields the schema leaves open are kept in',
            '//! `extra`. Regenerate with `python3 scripts/gen_kinds.py`.',
            '#![allow(non_snake_case, non_camel_case_types, clippy::all)]',
            '',
            'use std::collections::BTreeMap;',
            'use std::io::Write;',
            '',
            'use serde::{Deserialize, Serialize};',
            'use serde_json::Value;',
            '',
            'use crate::error::QvfError;',
            'use crate::reader::{Reader, Section, Tensor};',
            'use crate::typed::{',
            '    add_section, check_kind, json_member, opaque_member, optional, pattern_members, put_json, put_opaque, put_tensor,',
            '    section_fields, tensor_member, OpaqueData, TypedKind,',
            '};',
            'use crate::writer::{WriteError, Writer};',
            '']
    return '\n'.join(head + types.lines + items + any_kind), types.fallbacks
