const rye = require('rye');
const math = require('mathjs');
const gaussian = require('gaussian');
const sjcl = require('sjcl');
const tsUtils = require('./trapSamplerUtils.js'); 
const matUtils = require('./matUtils.js');
const PrimeField = rye.PrimeField;

class TrapSampler {
    constructor(q,n,m,sigma) {
        this.q = q;
        this.n = n;
        this.m = m;
        this.l = math.ceil(math.log(q,2));
        this.ZZq = new PrimeField(q);
        this.dist = new DiscreteGaussian(0,sigma);
        this.g = tsUtils.getGadgetVector(this.l)
        this.G = tsUtils.getGadget(this.g,this.n);
        this.zMap = tsUtils.createBuckets(this.q,this.l,this.g,this.dist,null);
        this.elen = this.m + this.n*this.l
        this.trapWidth = this.n*this.l
    }

    /**
     * Generates a uniform matrix AA along with a trapdoor matrix RR
     * @return {Array} [AA,RR]
     */
    genMatrixWithTrapdoor() {
        let A = math.mod(math.random([this.n,this.m]), this.q);
        let R = matUtils.genShortMatrix(this.m,this.trapWidth,this.dist);
        // console.log(math.multiply(A,R).valueOf());
        console.log(this.G.valueOf());
        let A1 = tsUtils.balance(math.mod(math.subtract(this.G, math.multiply(A,R)),this.q));
        let IR = math.eye(this.trapWidth);
        let AA = matUtils.augment(A,A1);
        let RR = matUtils.stack(R,IR);

        return [AA,RR]
    }

    /**
     * Return a pre-image x such that AA*X = u mod q for a vector u
     * @param {Matrix | Array} AA 
     * @param {Matrix | Array} RR 
     * @param {Vector} u 
     * @return {Vector} x
     */
    vecPreImage(AA,RR,u) {
        let pArr = []
        for (let i=0;i<this.elen;i++) {
            pArr.push(this.dist.sample);
        }   
        let p = matUtils.wrapVector(pArr);
        // v = (u-AA*p) mod q
        let v = (math.mod(math.subtract(u,math.multiply(AA,p)), this.q));
        try {
            let z = zPreImg(v);
        } catch(error) {
            console.error(error);
            return null;
        }
        // x = (p + RR*z) mod q 
        let x = math.mod(math.add(p,math.multiply(RR,z)),this.q);
        // res = (AA*x) mod q
        let res = math.mod(math.multiply(AA,x),this.q)
        if (!math.equal(res,u)) {
            throw new Error("Incorrect vector pre-image calculated, res: " + res + ", u: " + u);
        } else {
            return x;
        }
    }

    /**
     * For a given vector v, constructs a z such that 
     * G*z = v mod q
     * @param {Vector (mod q)} v
     * @return {Vector} z 
     */
    zPreImg(v) {
        let z = [];
        let vArr = v.valueOf();
        vArr.forEach(function(element) {
            let samps = this.zMap[element];
            if (matUtils.colSize(samps) != 0) {
                z.push(samps.pop());
            } else {
                console.log("There are no samples left for: ", + element);
                this.zMap = tsUtils.createBuckets(this.q, this.l, this.g, this.dist, this.zMap);
                throw new Error("New z values sampled, please retry.");
                // Not sure if this ever gets ran?
                return null;
            }
        });
        return matUtils.wrapVector(z);
    }
}

class DiscreteGaussian {
    constructor(mean,variance) {
        this.dist = gaussian(mean, variance);
    }
    get sample() {
        return math.round(this.dist.ppf(math.random()));
    }   
}

let m1 = math.matrix([[1,4,-3],[2,3,-4]]);
let m2 = math.matrix([[0,-2,-7],[5,0,0]]);
const ts = new TrapSampler(143,4,3,2);
// mats = ts.genMatrixWithTrapdoor();
// console.log(math.mod(math.multiply(mats[0],mats[1])),ts.q);
console.log(matUtils.tensorProduct(math.eye(4), m1));
