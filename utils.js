'no strict';
const math = require('mathjs');

const methods = {
    /**
     * Writes the value x mod q into the 
     * the bounds -q/2 < x < q/2
     * @param {Integer} x 
     * @param {Prime number} q 
     * @return {Number} xq
     */
    balance: function(x,q) {
        if (!q) {
            throw new Error("bad error");
        }
        let xq = math.mod(x,q);
        if (xq > math.floor(q/2)) {
            xq = xq - q;
        }
        return xq
    },
}
module.exports = methods;