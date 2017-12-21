/**
 * Example file for coela.js, also used for testing
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
const TrapSampler = require('./trapdoor.js');
const LWESampler = require('./lweSampler.js');
const GGH15Sampler = require('./ggh15Sampler.js');
const matUtils = require('./matUtils.js');
const math = require('mathjs');

class TrapdoorExample {
	constructor(q, n, m, sigma)	{
		this.ts = new TrapSampler(q,n,m,sigma);
	}

	generate() {
		let mats = this.ts.genMatrixWithTrapdoor();
		let U = math.matrix(math.round(math.multiply(math.random([this.ts.n, this.ts.n]),10)));
		let preImg = this.ts.matPreImage(mats.AA,mats.trap,U);

		return {mats: mats, U: U, preImg: preImg};
	}

	check(res) {
		let AA = res.mats.AA;
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
		let lwe = this.ls.shortSecretLWE(mats2.AA);
		let B = lwe.B;
		let preImg = this.ts.matPreImage(mats1.AA,mats1.trap,B);

		return {mats1: mats1, mats2: mats2, B: B, preImg: preImg};
	}

	check(res) {
		let AA0 = res.mats1.AA;
		let preImg = res.preImg;
		let B = res.B;
		let rowsB = matUtils.rowSize(B);
		let colsB = matUtils.colSize(B);

		return this.ts.checkMatrixPreImage(AA0,preImg,B,rowsB,colsB);
	}
}

class GGH15Example {
	constructor(q, n, m, sigma, kappa, branches, bookends) {
		this.ggh15 = new GGH15Sampler(q, n, m, sigma, kappa, branches, bookends);
		this.ls = this.ggh15.lwe;
		this.ts = this.ggh15.ts;
	}

	generate(width) {
		let nodes = this.ggh15.generateNodes();
		let encodings = this.ggh15.generateEdges(nodes.nodes, nodes.bookendMats, width);
		return {nodes: nodes, encodings: encodings};
	}

	check(res,branches,width) {

	}
}

let q = 143;
let n = 4;
let m = 3; 
let sigma = 2.0;

let trap = new TrapdoorExample(q, n, m, sigma);
let resTrap = trap.generate();
if (!trap.check(resTrap)) {
	throw new Error("Matrix pre-image check for Trapdoor example failed.");
}

let lwe = new LWEExample(q, n, m, sigma);
let resLWE = lwe.generate();
if (!lwe.check(resLWE)) {
	throw new Error("Matrix pre-image check for LWE example failed.");
}

let kappa = 3;
for (let branches=1;branches<3;branches++) {
	for (let width=1;width<3;width++) {
		let ggh15NoBookends = new GGH15Example(q, n, m, sigma, kappa, branches, false);
		let resGGH15NoBookends = ggh15NoBookends.generate(width);
		console.log(resGGH15NoBookends);

		let ggh15Bookends = new GGH15Example(q, n, m, sigma, kappa, branches, true);
		let resGGH15Bookends = ggh15Bookends.generate(width);
		console.log(resGGH15Bookends);
	}
}

