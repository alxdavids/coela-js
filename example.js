/**
 * Example file for coela.js, also used for testing
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
const TrapSampler = require('./trapdoor.js');
const LWESampler = require('./lweSampler.js');
const matUtils = require('./matUtils.js');
const math = require('mathjs');

class TrapdoorExample {
	constructor(q, n, m, sigma)	{
		this.ts = new TrapSampler(q,n,m,sigma);
	}

	generate() {
		let mats = this.ts.genMatrixWithTrapdoor();
		let U = math.matrix(math.round(math.multiply(math.random([this.ts.n, this.ts.n]),10)));
		let preImg = this.ts.matPreImage(mats[0],mats[1],U);

		return {mats: mats, U: U, preImg: preImg};
	}

	check(res) {
		let AA = res.mats[0];
		let preImg = res.preImg;
		let U = res.U;
		let rowsU = matUtils.rowSize(U);
		let colsU = matUtils.colSize(U);

		return this.ts.checkMatrixPreImage(AA,preImg,U,rowsU,colsU);
	}
}

class LWEExample {
	constructor(q, n, m, sigma) {
		this.ls = new LWESampler(q, n, m, sigma);
		this.ts = this.ls.ts;
	}

	generate() {
		let mats1 = this.ts.genMatrixWithTrapdoor();
		let mats2 = this.ts.genMatrixWithTrapdoor();
		let lwe = this.ls.sampleShortLWE(mats2[0]);
		let B = lwe.B;
		let preImg = this.ts.matPreImage(mats1[0],mats1[1],B);

		return {mats1: mats1, mats2: mats2, B: B, preImg: preImg};
	}

	check(res) {
		let AA0 = res.mats1[0];
		let preImg = res.preImg;
		let B = res.B;
		let rowsB = matUtils.rowSize(B);
		let colsB = matUtils.colSize(B);

		return this.ts.checkMatrixPreImage(AA0,preImg,B,rowsB,colsB);
	}
}

let trap = new TrapdoorExample(143, 4, 3, 2.0);
let resTrap = trap.generate();
if (!trap.check(resTrap)) {
	throw new Error("Matrix pre-image check for Trapdoor example failed.");
}

let lwe = new LWEExample(143, 4, 3, 2.0);
let resLWE = lwe.generate();
if (!lwe.check(resLWE)) {
	throw new Error("Matrix pre-image check for LWE example failed.");
}
