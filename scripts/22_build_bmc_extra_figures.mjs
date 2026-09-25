import fs from 'node:fs/promises';
import path from 'node:path';
import sharp from 'sharp';

const root=path.resolve(process.env.JOB49_ROOT || '.');
const out=path.join(root,'figures');
await fs.mkdir(out,{recursive:true});
const esc=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;');
const T=(x,y,v,size=28,weight=400,color='#243447',anchor='start')=>`<text x="${x}" y="${y}" fill="${color}" font-family="Arial, sans-serif" font-size="${size}" font-weight="${weight}" text-anchor="${anchor}">${esc(v)}</text>`;
const L=(x1,y1,x2,y2,color='#C6D2DC',width=2)=>`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="${width}"/>`;
const C=(x,y,r,fill,stroke='#FFFFFF')=>`<circle cx="${x}" cy="${y}" r="${r}" fill="${fill}" stroke="${stroke}" stroke-width="3"/>`;
const frame=(w,h,body)=>`<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"><rect width="${w}" height="${h}" fill="#FFFFFF"/>${body}</svg>`;
const parse=t=>t.trim().split(/\r?\n/).map(x=>x.split(','));

// Figure 2: GEO replication across independently defined marker panels.
const geo=parse(await fs.readFile(path.join(root,'manuscript','Supplementary_Table_S3_GEO_replication.csv'),'utf8'));
const rows=geo.slice(1).map(r=>({cohort:r[0],proxy:r[1],n:+r[2],rho:+r[3],p:+r[4]}));
const labels={myeloid_identity:'Myeloid identity',macrophage_state:'Macrophage state',tumor_glial_proxy:'Glial/tumor-like proxy'};
let b='';
b+=T(110,110,'Spearman correlation with the Hallmark xenobiotic score',30,400,'#475569');
const x0=750,x1=1910,lo=-0.5,hi=0.9,scale=v=>x0+(v-lo)/(hi-lo)*(x1-x0);
for(const tick of [-0.4,0,0.4,0.8]){let x=scale(tick);b+=L(x,250,x,1130,tick===0?'#64748B':'#DCE5EC',tick===0?3:2)+T(x,1195,tick.toFixed(1),25,400,'#475569','middle');}
const order=['GSE16011','GSE149009'];
let y=335;
for(const cohort of order){
  b+=T(105,y-65,`${cohort} (n=${rows.find(r=>r.cohort===cohort).n})`,30,700,'#183153');
  for(const key of ['myeloid_identity','macrophage_state','tumor_glial_proxy']){
    const r=rows.find(z=>z.cohort===cohort&&z.proxy===key);
    const color=key==='myeloid_identity'?'#2879A5':key==='macrophage_state'?'#38A38B':'#C57A50';
    b+=T(150,y,labels[key],27,400,'#243447');
    b+=L(scale(0),y,scale(r.rho),y,color,8)+C(scale(r.rho),y,15,color);
    b+=T(2030,y,`ρ = ${r.rho.toFixed(3)}`,28,700,color);
    y+=125;
  }
  y+=95;
}
b+=T(1330,1250,'Spearman ρ',28,400,'#475569','middle');
b+=T(110,1330,'The historical GEO labels do not establish IDH-wildtype eligibility.',25,400,'#657482');
await sharp(Buffer.from(frame(2400,1400,b))).png().withMetadata({density:300}).toFile(path.join(out,'Figure_2.png'));

// Figure 4: patient-paired single-cell contrast, preserving platform separation.
const sc=parse(await fs.readFile(path.join(root,'results','GSE131928_xenobiotic_patient_lineage.csv'),'utf8'));
const map=new Map();
for(const r of sc.slice(1)){
  const key=r[0]+'|'+r[1];if(!map.has(key))map.set(key,{});
  map.get(key)[r[2]]={cells:+r[3],score:+r[4]};
}
const groups={Smartseq2:[], '10X':[]};
for(const [key,v] of map){if(v.Myeloid&&v.Glial_tumor_like&&v.Myeloid.cells>=5&&v.Glial_tumor_like.cells>=5){const [platform,patient]=key.split('|');groups[platform].push({patient,diff:v.Myeloid.score-v.Glial_tumor_like.score});}}
for(const g of Object.values(groups))g.sort((a,b)=>a.diff-b.diff);
const test=parse(await fs.readFile(path.join(root,'results','GSE131928_xenobiotic_paired_tests.csv'),'utf8')).slice(1);
let q='';q+=T(110,110,'Myeloid minus glial/tumor-like mean pathway score',30,400,'#475569');
const panels=[{name:'Smartseq2',x:100,ys:285},{name:'10X',x:1270,ys:285}];
for(const panel of panels){
  const vals=groups[panel.name],xzero=panel.x+520,s=600;
  const stat=test.find(r=>r[0]===panel.name);
  q+=T(panel.x,panel.ys-50,`${panel.name} (paired patients=${vals.length})`,30,700,'#183153');
  q+=T(panel.x,panel.ys-5,`Wilcoxon P=${Number(stat[4]).toFixed(3)}`,25,400,'#64748B');
  q+=L(xzero,panel.ys+30,xzero,panel.ys+900,'#8394A4',3);
  for(let i=0;i<vals.length;i++){
    const y=panel.ys+95+i*83;const d=vals[i].diff,x=xzero+d*s;
    q+=L(panel.x+180,y,xzero+360,y,'#ECF0F4',2);
    q+=T(panel.x+10,y+8,vals[i].patient,25,400,'#344256');
    q+=L(xzero,y,x,y,d>=0?'#2879A5':'#C57A50',6)+C(x,y,12,d>=0?'#2879A5':'#C57A50');
    q+=T(panel.x+960,y+8,(d>=0?'+':'')+d.toFixed(3),24,400,'#475569','end');
  }
  q+=T(xzero,panel.ys+985,'0',25,400,'#475569','middle');
}
q+=T(120,1410,'Patient-paired differences vary in direction; the full program is not consistently myeloid enriched.',27,400,'#657482');
await sharp(Buffer.from(frame(2400,1500,q))).png().withMetadata({density:300}).toFile(path.join(out,'Figure_4.png'));
console.log(JSON.stringify({figure2:rows.length,figure4:groups.Smartseq2.length+groups['10X'].length}));
