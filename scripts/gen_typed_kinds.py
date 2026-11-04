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

NAMES={'wavefunction.gto':'WavefunctionGTO','wavefunction.bloch':'WavefunctionBloch','job.spec':'JobSpec','run.record':'RunRecord','dos.total':'DosTotal','dos.projected':'DosProjected','dos.coop':'DosCoop','dos.cohp':'DosCohp','topology.qtaim':'TopologyQTAIM','spectra.ir':'SpectraIR','spectra.raman':'SpectraRaman','spectra.uvvis':'SpectraUVVis','spectra.ecd':'SpectraECD','spectra.vcd':'SpectraVCD','spectra.nmr':'SpectraNMR','spectra.epr':'SpectraEPR','spectra.xas':'SpectraXAS','spectra.xes':'SpectraXES','spectra.mossbauer':'SpectraMossbauer','basis.ao':'BasisAO','volume.elf':'VolumeELF','volume.rdg':'VolumeRDG','scf_history':'ScfHistory'}
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
