/**
 * Parameter selection (primarily for GGH15 purposes)
 *
 */

/*jshint esversion: 6 */
/*jshint node: true */
const math = require('mathjs');
const sprintf = require('sprintf-js').sprintf;
const methods = {
    /**
     * Choose parameters according to Appendix A of GGH15
     * @param {BigNumber} lambda security parameter 
     * @param {BigNumber} n lwe dimension
     * @param {BigNumber} q modulus
     * @param {BigNumber} kappa multilinearity
     * @returns {Literal}
     */
    getParameters(lambda,kappa) {
        // helpers
        let kl = math.multiply(kappa, lambda);
        let big2 = math.bignumber(2);
        let logkl = math.log(kl, big2);
        let k_1 = math.subtract(kappa,math.bignumber(1));
        
        // LWE parameter selection
        let n = math.floor(math.pow(lambda,lambda));
        let q = math.floor(math.pow(big2, math.divide(n, lambda)));
        console.log(sprintf("n:%v q:%v",n,q));
        const s = math.sqrt(n);
        let chkq = math.floor(math.multiply(s, math.pow(big2, math.divide(n, lambda))));
        if (math.compare(q, chkq) == 1) {
            throw new Error(sprintf("Parameter selection is incorrect for q, required: < %v, given: %v", chkq, q));
        }

        // trapdoor parameter selection
        const m = math.multiply(n,math.log(q,big2));
        const sigma = math.sqrt(math.multiply(n,math.log(q,big2)));

        // check kappa selection
        const coeff = math.multiply(math.pow(n,kappa),math.pow(m,k_1));
        const qpow = math.multiply(kappa,math.pow(math.log(q,big2),math.divide(k_1,big2)));
        console.log(sprintf("coeff:%v qpow:%v",coeff,qpow));
        let chkq2 = math.floor(math.multiply(coeff,qpow));
        let q34 = math.floor(math.pow(q, math.divide(math.bignumber(3), math.bignumber(4))));
        if (math.compare(chkq2,q34) == 1) {
            throw new Error(sprintf("Parameter selection is incorrect for kappa, required: n^kappa x m^(kappa-1) x kappa x (log q)^((kappa-1)/2) < q^(3/4), LHS: %v, RHS: %v", chkq2, q34));
        }

        return {q: q, n: n, s: s, m: m, sigma: sigma, lambda: lambda, kappa: kappa};
    }
}
module.exports = methods.getParameters;
