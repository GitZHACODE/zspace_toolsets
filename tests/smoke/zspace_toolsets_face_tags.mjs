// Node-only integration check against the locally built Toolsets runtime.
// Usage: node tests/smoke/zspace_toolsets_face_tags.mjs <input.usda> [mode=0]
import assert from 'node:assert/strict';
import {readFile,writeFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import factory from '../../wasm/out/zspace_toolset.js';
const module=await factory({slicingThreadCount:4,wasmBinary:await readFile(new URL('../../wasm/out/zspace_toolset.wasm',import.meta.url)),print(){}});
const call=(name,args=[])=>module.ccall(name,'number',args.map(v=>typeof v==='string'?'string':'number'),args);
const check=n=>{assert(n,module.UTF8ToString(call('zspace_toolsets_last_error')));return n;};
const upload=(h,slot,a)=>{
  const buffer=Float64Array.from(a),p=module._zspace_alloc(buffer.byteLength);
  try {new Float64Array(module.HEAPU32.buffer,p,buffer.length).set(buffer);check(call('zspace_slicing_set_array',[h,slot,p,buffer.length]));}
  finally{module._zspace_free(p);}
};
const read=(h,r,s)=>{
 const n=call('zspace_slicing_array_count',[h,r,s]);assert(n>=0);
 const p=call('zspace_slicing_array_ptr',[h,r,s]);return Array.from(new Float64Array(module.HEAPU32.buffer,p,n));
};
let h;
try{
 const bytes=await readFile(process.argv[2]);module.FS.writeFile('/input.usda',bytes);
 h=check(call('zspace_slicing_create'));check(call('zspace_slicing_import_usd',[h,'/input.usda']));
 const input=Array.from({length:11},(_,slot)=>{
   const n=call('zspace_slicing_input_count',[h,slot]),p=call('zspace_slicing_input_ptr',[h,slot]);
   return Array.from(new Float64Array(module.HEAPU32.buffer,p,n));
 });
 const options=input[8];
 const fixture=[...input.slice(0,8),input[9],input[10],
   [0,1,4,2,6,7,5,3,13,14,15,16,10,11].map(i=>options[i])];
 await writeFile('build/shared-facetags-input.txt',fixture.map(a=>a.length+'\n'+a.join(' ')).join('\n')+'\n');
 check(call('zspace_slicing_compute',[h,Number(process.argv[3]??0),1,1]));
 const paths=[];let combined;
 for(let r=0;r<call('zspace_slicing_record_count',[h]);++r){
  const kind=call('zspace_slicing_record_kind',[h,r]);if(kind!==1&&kind!==4)continue;
  const a=Array.from({length:7},(_,slot)=>read(h,r,slot));
  assert.equal(a[6].length,a[0].length/3,'FaceTag output lost alignment');
  assert(a[6].every(t=>t===0||[3,4,5,6].includes(t)));
  const path={layer:call('zspace_slicing_record_layer',[h,r]),a};
  if(kind===1)combined=path;else paths.push(path);
 }
 const counts=path=>Object.fromEntries([3,4,5,6].map(t=>[t,path.a[6].filter(v=>v===t).length]));
 const expected=counts(paths[0]);for(const p of paths)assert.deepEqual(counts(p),expected,'Per-point FaceTag counts changed across layers');
 assert(Object.values(expected).every(v=>v>0));
 check(call('zspace_slicing_export_begin',[h,paths.length]));
 for(const p of [...paths,combined]){
  for(const slot of [0,1,2,3,4,6])upload(h,slot,p.a[slot]);
  check(call('zspace_slicing_export_add_path',[h,p.layer]));
 }
 check(call('zspace_slicing_export_usd',[h,'/input.usda','/output.usda',0]));
 const output=module.FS.readFile('/output.usda',{encoding:'utf8'});
 assert(output.includes('printFaceTags'),'USD export omitted per-point FaceTags');
 const exported=output.match(/int\[\]\s+printFaceTags\s*=\s*\[([^\]]*)\]/);
 assert(exported,'Missing exported integer FaceTag array');
 assert.deepEqual(exported[1].split(',').map(v=>Number(v.trim())),combined.a[6],'USD export reordered FaceTags');
 await writeFile('build/facetags-output.usda',output);
 // Actual computed geometry for the Alice2 documentation sketch; no rendering here.
 const rows=[paths.length];
 for(const p of paths){rows.push(p.layer,p.a[0].length/3,p.a[1].length);for(let i=0;i<p.a[6].length;++i)rows.push(...p.a[0].slice(3*i,3*i+3),p.a[6][i]);rows.push(...p.a[1]);}
 if(Number(process.argv[3]??0)===0)await writeFile('build/facetags-review.dat',rows.join(' ')+'\n');
 const summary={input:process.argv[2],sha256:createHash('sha256').update(bytes).digest('hex'),mode:Number(process.argv[3]??0),layers:paths.length,pointCounts:expected,exportedFaceTags:true};
 await writeFile(`build/facetags-wasm-verification-${summary.mode}.json`,JSON.stringify(summary,null,2)+'\n');console.log(JSON.stringify(summary));
}finally{if(h)call('zspace_slicing_destroy',[h]);module.PThread.terminateAllThreads();}
