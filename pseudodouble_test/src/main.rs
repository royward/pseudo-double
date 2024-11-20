mod pseudo_double;

use pseudo_double::PseudoDouble;
use libm::{ ldexp };
use rand::{Rng,SeedableRng};
use rand::rngs::StdRng;

pub const NEAR_EXACT14: f64=0.99999999999999;
pub const NEAR_EXACT13: f64=0.9999999999999;
pub const NEAR_EXACT12: f64=0.999999999999;
pub const NEAR_EXACT11: f64=0.99999999999;
pub const NEAR_EXACT10: f64=0.9999999999;
pub const NEAR_EXACT9:  f64=0.999999999;
pub const NEAR_EXACT8:  f64=0.99999999;
pub const NEAR_EXACT7:  f64=0.9999999;
pub const NEAR_EXACT6:  f64=0.999999;
pub const NEAR_EXACT5:  f64=0.99999;
pub const NEAR_EXACT4:  f64=0.9999;
pub const NEAR_EXACT3:  f64=0.999;

const fn compare(d1:f64, d2:f64, exactness:f64) -> bool {
	if d1>=0.0 {
		return d1*exactness<=d2 && d2*exactness<=d1;
	} else {
		return d1*exactness>=d2 && d2*exactness>=d1;
	}
}

fn main() {
	let mut count=0u32;
	let mut failures=0u32;
	let mut list = Vec::new();
	for i in -20..20 {
		list.push(ldexp(1.0f64,i));
		list.push(-ldexp(1.0,i));
		list.push(ldexp(3.0,i));
		list.push(-ldexp(3.0,i));
		list.push(f64::from(i));
		list.push(f64::from(i)+0.5);
	}
	let mut rng = StdRng::seed_from_u64(222); // rand::thread_rng();
	for _i in 0..100 {
		let r:f64=rng.gen();
		let f1=r*1000000.0;
		list.push(f1);
		list.push(-f1);
		let f2=r/1000000.0;
		list.push(f2);
		list.push(-f2);
	}
	// {
	// 	let f1=9.5367431640625e-07;
	// 	let f2=5.4088650980917826e-07;
	// 	let pd1=PseudoDouble::double_to_pseudodouble(f1);
	// 	let pd2=PseudoDouble::double_to_pseudodouble(f2);
	// 	let b=pd1<pd2;
	// 	println!("compccc {}<{} {}<{}",f1,f2,f1<f2,b);
	// }
	for f1v in &list {
		let f1=*f1v;
		let pd1=PseudoDouble::double_to_pseudodouble(f1);
		for f2v in &list {
			let f2=*f2v;
			let pd2=PseudoDouble::double_to_pseudodouble(f2);
			let ffa=f64::from(pd1+pd2);
			count+=1;
			if !compare(f1+f2,ffa,NEAR_EXACT8) {
				failures+=1;
				println!("add {}+{}=={}!={}",f64::from(f1),f64::from(f2),f64::from(f1+f2),f64::from(ffa));
			}
			let ffs=f64::from(pd1-pd2);
			count+=1;
			if !compare(f1-f2,ffs,NEAR_EXACT8) {
				failures+=1;
				println!("sub {}-{}=={}!={}",f64::from(f1),f64::from(f2),f64::from(f1-f2),f64::from(ffs));
			}
			let ffm=f64::from(pd1*pd2);
			count+=1;
			if !compare(f1*f2,ffm,NEAR_EXACT8) {
				failures+=1;
				println!("mul {}*{}=={}!={}",f64::from(f1),f64::from(f2),f64::from(f1*f2),f64::from(ffm));
			}
			count+=1;
			if(PseudoDouble::max(pd1,pd2)==pd1) != (f64::max(f1,f2)==f1) {
				failures+=1;
				println!("difference in max({},{}) {} {}",f1,f2,f64::from(PseudoDouble::max(pd1,pd2)),f64::max(f1,f2));
				let t=PseudoDouble::max(pd1,pd2);
				println!("{}",f64::from(t));
			}
			count+=1;
			if(PseudoDouble::min(pd1,pd2)==pd1) != (f64::min(f1,f2)==f1) {
				failures+=1;
				println!("difference in min({},{})",f1,f2);
			}
			if f2!=0.0 {
				let ffd=f64::from(pd1/pd2);
				count+=1;
				if !compare(f1/f2,ffd,NEAR_EXACT8) {
					failures+=1;
					println!("div {}/{}=={}!={}",f64::from(f1),f64::from(f2),f64::from(f1/f2),f64::from(ffd));
				}
			}
			count+=1;
			if(f1<f2) != (pd1<pd2) {
				failures+=1;
				println!("comp {}<{} {}<{}",f1,f2,f1<f2,pd1<pd2);
				println!("{}",pd1<pd2);
			}
			count+=1;
			// println!("here atan2({},{})",f1,f2);
			let ft=pd1.atan2(pd2);
			if( !compare(f1.atan2(f2),f64::from(ft),NEAR_EXACT8)) && ((f2/f1).abs()<1e9 || !compare(f1.atan2(f2),f64::from(ft),NEAR_EXACT3)) && (f1>=0.0 || f2!=0.0) {
				failures+=1;
				println!("atan2({},{}) {} {}",f1,f2,f1.atan2(f2),f64::from(ft));
			}
		}
	}
	for i in -20..20 {
		list.push(f64::from(i)+0.4999999);
		list.push(f64::from(i)+0.5000001);
	}
	for fv in &list {
		let f=*fv;
		let pd=PseudoDouble::double_to_pseudodouble(f);
		let ffc=f64::from(pd);
		count+=1;
		if !compare(f,ffc,NEAR_EXACT13) {
			failures+=1;
			println!("conv {} {}",f,ffc);
		}
		let ffn=f64::from(-pd);
		count+=1;
		if !compare(-f,ffn,NEAR_EXACT13) {
			failures+=1;
			println!("neg {} {} {}",f,-f,ffn);
			let t=-pd;
			let tt=f64::from(t);
			println!("{}",tt)
		}
		if f>0.0 {
			count+=1;
			let ff=pd.inv_sqrt();
			if !compare(1.0/f.sqrt(),f64::from(ff),NEAR_EXACT13) {
				failures+=1;
				println!("inv_sqrt {} {}",1.0/f.sqrt(),f64::from(ff));
			}
		}
		if f>0.0 {
			count+=1;
			let ff=pd.sqrt();
			if !compare(f.sqrt(),f64::from(ff),NEAR_EXACT13) {
				failures+=1;
				println!("sqrt {} {}",f.sqrt(),f64::from(ff));
			}
		}
		let ff=pd.floor();
		count+=1;
		if !compare(f.floor(),f64::from(ff),NEAR_EXACT13) {
			failures+=1;
			println!("floor {} {} {}",f,f.floor(),f64::from(ff));
		}
		let fc=pd.ceil();
		count+=1;
		if !compare(f.ceil(),f64::from(fc),NEAR_EXACT13) {
			failures+=1;
			println!("ceil {} {} {}",f,f.ceil(),f64::from(fc));
		}
		let fr=pd.round();
		count+=1;
		if !compare(f.round(),f64::from(fr),NEAR_EXACT13) {
			failures+=1;
			println!("round {} {} {}",f,f.round(),f64::from(fr));
		}
		if f<128.0 && f>-128.0 {
			count+=1;
			let ff=pd.exp2();
			if !compare(f.exp2(),f64::from(ff),NEAR_EXACT12) {
				failures+=1;
				println!("exp2 {} {} {}",f,f.exp2(),f64::from(ff));
			}
		}
		if f<96.0 && f>-96.0 {
			count+=1;
			let ff=pd.exp();
			if !compare(f.exp(),f64::from(ff),NEAR_EXACT12) {
				failures+=1;
				println!("exp {} {} {}",f,f.exp(),f64::from(ff));
			}
		}
		if f>0.0 {
			count+=1;
			let ff=pd.log2();
			if !compare(f.log2(),f64::from(ff),NEAR_EXACT12) {
				failures+=1;
				println!("log2 {} {} {}",f,f.log2(),f64::from(ff));
			}
		}
		if -10000.0<f && f<10000.0 {
			let fs=pd.sin();
			count+=1;
			if !compare(f.sin(),f64::from(fs),NEAR_EXACT9) {
				failures+=1;
				println!("sin {} {} {}",f,f.sin(),f64::from(fs));
			}
			let fc=pd.cos();
			count+=1;
			if !compare(f.cos(),f64::from(fc),NEAR_EXACT9) {
				failures+=1;
				println!("cos {} {} {}",f,f.cos(),f64::from(fc));
			}
		}
	}
	for i in -1000i64..1000i64 {
		let pd=PseudoDouble::from(i);
		let ii=i64::from(pd);
		count+=1;
		if i!=ii {
			failures+=1;
			println!("i64 convert {} {}",i,ii);
		}
	}
	for i in 0u64..1000u64 {
		let pd=PseudoDouble::from(i);
		let ii=u64::from(pd);
		count+=1;
		if i!=ii {
			failures+=1;
			println!("u64 convert {} {}",i,ii);
		}
	}
	println!("Tests done, passed {}/{}",count-failures,count);

	let x=PseudoDouble::from(2);
	let y=PseudoDouble::from(7);
	let z=i64::from(x+y);
	println!("z={}",z);
}

