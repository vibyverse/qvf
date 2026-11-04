/**
 * Encoding of binary members for the writer (spec § 2.2, § 4.2): numbers,
 * nested arrays and typed arrays to flat, C-ordered, little-endian bytes of a
 * QVF dtype.
 */

import { isDType } from "./binary.js";
import { QvfError } from "./errors.js";
import { DTYPE_ITEMSIZE, type DType, type TypedArray } from "./model.js";

const HOST_IS_LITTLE_ENDIAN = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;

/**
 * Values a binary member can be written from: a typed array (flat, C-ordered),
 * or a possibly nested array of numbers or bigints whose nesting gives the
 * shape (it must be rectangular).
 */
export type TensorInput = TypedArray | NestedNumbers;

/** A number, a bigint, or an arbitrarily nested array of them. */
export type NestedNumbers = number | bigint | readonly NestedNumbers[];

/** The dtype a typed array stores, or `undefined` for a plain array. */
export function dtypeOf(data: unknown): DType | undefined {
  if (data instanceof Float64Array) return "float64";
  if (data instanceof Float32Array) return "float32";
  if (data instanceof Int8Array) return "int8";
  if (data instanceof Int16Array) return "int16";
  if (data instanceof Int32Array) return "int32";
  if (data instanceof BigInt64Array) return "int64";
  if (data instanceof Uint8Array) return "uint8";
  if (data instanceof Uint16Array) return "uint16";
  if (data instanceof Uint32Array) return "uint32";
  if (data instanceof BigUint64Array) return "uint64";
  return undefined;
}

function isTypedArray(data: unknown): data is TypedArray {
  return dtypeOf(data) !== undefined;
}

/**
 * Flatten a rectangular nested array into its values and shape. A scalar has
 * shape `[]`; a typed array inside a nested array is one innermost row.
 */
export function flatten(data: NestedNumbers | TypedArray, what = "data"): {
  values: ArrayLike<number | bigint>[];
  shape: number[];
} {
  const shape: number[] = [];
  let probe: unknown = data;
  while (Array.isArray(probe)) {
    shape.push(probe.length);
    if (probe.length === 0) break;
    probe = probe[0];
  }
  if (isTypedArray(probe)) shape.push(probe.length);
  const chunks: ArrayLike<number | bigint>[] = [];
  const walk = (value: unknown, depth: number, where: string): void => {
    if (depth === shape.length) {
      if (typeof value !== "number" && typeof value !== "bigint") {
        throw new TypeError(`${what}${where}: expected a number, got ${describe(value)}`);
      }
      chunks.push([value]);
      return;
    }
    if (isTypedArray(value) && depth === shape.length - 1) {
      if (value.length !== shape[depth]) {
        throw new RangeError(`${what}${where}: ragged array, length ${value.length} where ${shape[depth]} is expected`);
      }
      chunks.push(value as ArrayLike<number | bigint>);
      return;
    }
    if (!Array.isArray(value)) {
      throw new TypeError(`${what}${where}: expected an array, got ${describe(value)}`);
    }
    if (value.length !== shape[depth]) {
      throw new RangeError(`${what}${where}: ragged array, length ${value.length} where ${shape[depth]} is expected`);
    }
    if (depth === shape.length - 1) {
      for (let i = 0; i < value.length; i++) {
        const v: unknown = value[i];
        if (typeof v !== "number" && typeof v !== "bigint") {
          throw new TypeError(`${what}${where}[${i}]: expected a number, got ${describe(v)}`);
        }
      }
      chunks.push(value as (number | bigint)[]);
      return;
    }
    for (let i = 0; i < value.length; i++) walk(value[i], depth + 1, `${where}[${i}]`);
  };
  walk(data, 0, "");
  return { values: chunks, shape };
}

function describe(value: unknown): string {
  if (value === null) return "null";
  if (Array.isArray(value)) return "an array";
  return typeof value;
}

function product(shape: readonly number[]): number {
  let n = 1;
  for (const d of shape) n *= d;
  return n;
}

const INTEGER_RANGE: Partial<Record<DType, [number, number]>> = {
  int8: [-128, 127],
  int16: [-32768, 32767],
  int32: [-2147483648, 2147483647],
  uint8: [0, 255],
  uint16: [0, 65535],
  uint32: [0, 4294967295],
};

const BIG_RANGE: Partial<Record<DType, [bigint, bigint]>> = {
  int64: [-(1n << 63n), (1n << 63n) - 1n],
  uint64: [0n, (1n << 64n) - 1n],
};

/** Write one value at element `i` of `view` (little-endian). Throws for a value the dtype cannot hold. */
function put(view: DataView, dtype: DType, i: number, value: number | bigint, what: string): void {
  const o = i * DTYPE_ITEMSIZE[dtype];
  if (dtype === "float32" || dtype === "float64") {
    const x = typeof value === "bigint" ? Number(value) : value;
    if (dtype === "float32") view.setFloat32(o, x, true);
    else view.setFloat64(o, x, true);
    return;
  }
  const big = BIG_RANGE[dtype];
  if (big) {
    if (typeof value === "number" && !Number.isSafeInteger(value)) {
      throw new RangeError(`${what}: element ${i} (${value}) is not a safe integer for ${dtype}; pass a bigint`);
    }
    const b = BigInt(value);
    if (b < big[0] || b > big[1]) throw new RangeError(`${what}: element ${i} (${b}) is out of range for ${dtype}`);
    if (dtype === "int64") view.setBigInt64(o, b, true);
    else view.setBigUint64(o, b, true);
    return;
  }
  const x = typeof value === "bigint" ? Number(value) : value;
  const range = INTEGER_RANGE[dtype]!;
  if (!Number.isInteger(x) || x < range[0] || x > range[1]) {
    throw new RangeError(`${what}: element ${i} (${x}) is not an integer in the ${dtype} range`);
  }
  switch (dtype) {
    case "int8":
      view.setInt8(o, x);
      break;
    case "uint8":
      view.setUint8(o, x);
      break;
    case "int16":
      view.setInt16(o, x, true);
      break;
    case "uint16":
      view.setUint16(o, x, true);
      break;
    case "int32":
      view.setInt32(o, x, true);
      break;
    case "uint32":
      view.setUint32(o, x, true);
      break;
  }
}

/** An encoded binary member: its little-endian bytes, dtype and shape. */
export interface EncodedTensor {
  bytes: Uint8Array;
  dtype: DType;
  shape: number[];
}

/**
 * Encode `data` as the little-endian bytes of `dtype` (spec § 2.2).
 *
 * `dtype` defaults to the typed array's own dtype, and to `float64` for plain
 * arrays. `shape` defaults to the nesting of a plain array, or `[length]` for
 * a typed array; a given shape must hold exactly the number of elements.
 * A typed array of the requested dtype on a little-endian host is copied
 * once; anything else is converted element by element, and an integer dtype
 * refuses values it cannot hold exactly.
 */
export function encodeTensor(
  data: TensorInput,
  options: { dtype?: DType; shape?: readonly number[]; what?: string } = {},
): EncodedTensor {
  const what = options.what ?? "data";
  const own = dtypeOf(data);
  const dtype = options.dtype ?? own ?? "float64";
  if (!isDType(dtype)) throw new QvfError("E_DTYPE", `${what}: ${JSON.stringify(dtype)} is not a QVF dtype (spec § 4.2)`);
  const itemsize = DTYPE_ITEMSIZE[dtype];
  let flat: ArrayLike<number | bigint>[];
  let inferred: number[];
  if (own !== undefined) {
    flat = [data as ArrayLike<number | bigint>];
    inferred = [(data as TypedArray).length];
  } else {
    const f = flatten(data as NestedNumbers, what);
    flat = f.values;
    inferred = f.shape;
  }
  const count = flat.reduce((n, c) => n + c.length, 0);
  let shape = inferred;
  if (options.shape !== undefined) {
    shape = [...options.shape];
    if (!shape.every((d) => Number.isInteger(d) && d >= 0)) {
      throw new QvfError("E_MEMBER_MALFORMED", `${what}: shape ${JSON.stringify(shape)} must be non-negative integers`);
    }
    if (product(shape) !== count) {
      throw new QvfError(
        "E_BINARY_SIZE",
        `${what}: shape ${JSON.stringify(shape)} holds ${product(shape)} elements, the data has ${count}`,
      );
    }
  }
  const bytes = new Uint8Array(count * itemsize);
  if (own === dtype && HOST_IS_LITTLE_ENDIAN) {
    const t = data as TypedArray;
    bytes.set(new Uint8Array(t.buffer, t.byteOffset, t.byteLength));
    return { bytes, dtype, shape };
  }
  const view = new DataView(bytes.buffer);
  let i = 0;
  for (const chunk of flat) {
    for (let k = 0; k < chunk.length; k++) put(view, dtype, i++, chunk[k]!, what);
  }
  return { bytes, dtype, shape };
}
