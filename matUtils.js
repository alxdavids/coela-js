/*
 * Matrix utility functions
 */
const math = require('mathjs');

const methods = {
    /**
     * Get the col size of the matrix
     * @param {Matrix | Array} matrix
     * @return {Integer} size
     */
    rowSize: function(matrix) {
        let matrixSize = math.size(matrix).valueOf();
        if (!matrixSize[1]) {
            return 1;
        } else {
            return math.size(matrix).valueOf()[0];
        }
    },

    /**
     * Get the row size of the matrix
     * @param {Matrix | Array} matrix
     * @return {Integer} size
     */
    colSize: function(matrix) {
        let matrixSize = math.size(matrix).valueOf();
        if (!matrixSize[1]) {
            return matrixSize[0];
        } else {
            return math.size(matrix).valueOf()[1];
        }
    },

    /**
     * Wrap an array into a 1D vector
     * @param {Array} arr 
     * @return {Matrix | Array} Returns the array as a 1D vector
     */
    wrapVector: function(arr) {
        mat = [];
        mat.push(arr);
        return math.matrix(arr);
    },

    /**
    * Retrieve a column from a matrix
    * @param {Matrix | Array} matrix
    * @param {number} index    Zero based column index
    * @return {Matrix | Array} Returns the column as a vector
    */
    getCol: function(matrix, index) {
        let rows = rowSize(matrix);
        if (index == 0 && methods.colSize(matrix) == 1) {
            return math.flatten(matrix.valueOf());
        } else if (methods.colSize(matrix) == 1) {
            throw new Error("Index is non-zero but number of cols is 1");
        } else {
            return math.flatten(math.subset(matrix, math.index(math.range(0,rows), index)));
        }
    },

    /**
    * Retrieve a row from a matrix
    * @param {Matrix | Array} matrix
    * @param {number} index    Zero based row index
    * @return {Matrix | Array} Returns the row as a vector
    */
    getRow: function(matrix, index) {
        let cols = methods.colSize(matrix);
        if (index == 0 && methods.rowSize(matrix) == 1) {
            return math.flatten(matrix.valueOf());
        } else if (methods.rowSize(matrix) == 1) {
            throw new Error("Index is non-zero but number of rows is 1");
        } else {
            return math.flatten(math.subset(matrix, math.index(index, math.range(0,cols))));
        }
    },

    /**
     * Compute the tensor product of two matrices
     * @param {Matrix | Array} m1 
     * @param {Matrix | Array} m2 
     * @return {Matrix | Array} Returns the resulting matrix
     */
    tensorProduct: function(m1,m2) {
        let mTensor = [];
        let m2RowSize = methods.rowSize(m2);
        let m1Arr = m1.valueOf();
        m1Arr.forEach(function(row,index){
            let matRow = [];
            row.forEach(function(entry,index) {
                matRow.push(math.multiply(entry,m2));
            });
            for (let k=0; k<m2RowSize; k++) {
                let row = [];
                matRow.forEach(function(matrix,index) {
                    row = row.concat(methods.getRow(matrix,k).valueOf());
                });
                mTensor.push(row);
            }
        });
        return math.matrix(mTensor);
    },

    /**
     * Concatenate a matrix m2 to the RHS of m1
     * @param {Matrix | Array} m1 
     * @param {Matrix | Array} m2 
     * @return {Matrix | Array} M = [m1|m2]
     */
    augment: function(m1,m2) {
        let newRows = [];
        m1.valueOf().forEach(function(row,index) {
            row = row.concat(methods.getRow(m2,index).valueOf());
            newRows.push(row);
        });
        return math.matrix(newRows);
    },

    /**
     * Concatenate a matrix m2 below m1
     * @param {Matrix | Array} m1 
     * @param {Matrix | Array} m2 
     * @return {Matrix | Array} M = [m1]
     *                              ----
     *                              [m2]
     */
    stack: function(m1,m2) {
        let m1Arr = m1.valueOf();
        let m2Arr = m2.valueOf();
        m2Arr.forEach(function(row,index) {
            m1Arr.push(row);
        });
        return math.matrix(m1Arr);
    },

    /**
     * Generates a matrix with Gaussian-sampled entries
     * @param {Number} rows 
     * @param {Number} cols 
     * @return {Matrix | Array} R
     */
    genShortMatrix: function(rows,cols,dist) {
        matrix = [];
        for (let i=0; i<rows; i++) {
            row = [];
            for (let j=0; j<cols; j++) {
                row.push(dist.sample);
            }
            matrix.push(row);
        }
        return math.matrix(matrix);
    }
}
module.exports = methods;