#ifndef BINMAT_H
#define BINMAT_H
//#include <ostream>
#include <iostream>
#include <bitset>
#include <cstring>

typedef unsigned long idx_t;
typedef unsigned long block_t;
//typedef unsigned char block_t; // for debugging!
//typedef unsigned short block_t; // for debugging!

#define BITS_PER_BLOCK (sizeof(block_t)*8)
#define ONES   (~block_t(0))
#define ZEROES (block_t(0))
#define LSB    block_t(1)
#define MSB    (LSB<<(BITS_PER_BLOCK-1))
#define IMSB   (ONES>>1)
#define ILSB   (ONES<<1)

#define XOR(a,b) ( (!(a) && (b)) ||  ((a) && !(b)) )

typedef std::bitset<BITS_PER_BLOCK> bm_bitset; // for printing

/** 
 * ROW major (C-style) packed binary matrix 
 * Whenever possible, operations are carried out at block level.
 */
class binary_matrix {
public:

  /**
   * Allocates a binary matrix of the specified dimensions.
   */
  binary_matrix(idx_t _rows, idx_t _cols);

 binary_matrix() { memset(this,0,sizeof(binary_matrix)); }

 /**
  * Deep copy, with allocation
  */
 binary_matrix(const binary_matrix& other);

  /** 
   * Destructor
   */
  ~binary_matrix() {  } // we delete explicitly

 void allocate(idx_t _rows, idx_t _cols);

 binary_matrix get_vectorized() const;

 binary_matrix get_col(const idx_t j) const;

 binary_matrix get_row(const idx_t i) const ;

 binary_matrix get_submatrix(const idx_t i0,const  idx_t i1,const  idx_t j0,const  idx_t j1) const;

 binary_matrix get_copy() const;

 binary_matrix get_transposed() const;

 void copy_vectorized_to(binary_matrix& B) const;

 void copy_col_to(const idx_t j, binary_matrix& B) const;

 void copy_row_to(const idx_t i, binary_matrix& B) const ;

 void copy_submatrix_to(const idx_t i0,const  idx_t i1,const  idx_t j0,const  idx_t j1, binary_matrix& B) const;

 void copy_to(binary_matrix& B) const;

 void transpose_to(binary_matrix& B) const;

 void set_vectorized(const binary_matrix& src);

 void set_col(const idx_t j, const binary_matrix& src);

 void set_row(const idx_t i, const binary_matrix& src);

 void set_submatrix(const idx_t i0,const  idx_t j0, const binary_matrix& src);

 void add_rows(idx_t nrows);

 void remove_rows(idx_t nrows); // lightweight, does not deallocate space!

 /** @return the number of rows of the matrix */
 inline idx_t get_rows() const { return rows; }

 /** @return the number of columns of the matrix */
 inline idx_t get_cols() const { return cols; }

 /** @return the total number of bits of the matrix */
 inline idx_t get_len() const { return len; }

 /**
  * Sets all bits to zero. FIX: will erase whole block!
  */
 void clear(); 

 /**
  * Sets all bits to one.
  */
 void set();

 /**
  * Flip bits
  */
 void flip();

 /**
  * @return the bit at the specified position
  */
 inline bool get(const idx_t i, const idx_t j) const { 
   return data[i*blocks_per_row + (j / BITS_PER_BLOCK)] & (MSB >> (j % BITS_PER_BLOCK)); 
 }

 /**
  * Sets the bit at the specified position
  */
 inline void set(const idx_t i, const idx_t j, const bool v) {
   if (v) set(i,j); else clear(i,j);
 }

 /**
  * Sets the bit at the specified position to 1
  */
 inline void set(const idx_t i, const idx_t j) {
   data[i*blocks_per_row + (j / BITS_PER_BLOCK)] |= (MSB >> (j % BITS_PER_BLOCK));
 }

 inline void flip(const idx_t i, const idx_t j) {
   data[i*blocks_per_row + (j / BITS_PER_BLOCK)] ^= (MSB >> (j % BITS_PER_BLOCK));
 }

 /**
  * Sets the bit at the specified position to 0
  */
 inline void clear(const idx_t i, const idx_t j) {
   data[i*blocks_per_row + (j / BITS_PER_BLOCK)] &= ONES ^ (MSB >> (j % BITS_PER_BLOCK)); 
 }

 /** @return the number of ones in the matrix */
 idx_t weight() const; 

 /** @return the number of ones in the given row */
 idx_t row_weight(idx_t i) const; 

 /** @return the number of ones in the given column */
 idx_t col_weight(idx_t j) const; 

 /** @return the boolean sum (XOR), that is, the parity, of the whole matrix */
 bool sum() const; 

 /** @return the boolean sum (XOR), that is, the parity, of a given row */
 bool row_sum(idx_t i) const; 

 /** @return the boolean sum (XOR), that is, the parity, of a given column */
 bool col_sum(idx_t j) const; 

 /** Overloaded print operator */
 friend std::ostream& operator<<(std::ostream& out, const binary_matrix& A);

 /** C = A XOR B. It is assumed that C contains enough space for the result. */
 friend binary_matrix& add(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

#define bool_xor add

 /** C = A AND B. It is assumed that C contains enough space for the result. */
 friend binary_matrix& bool_and(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** 
  * Linear matrix product, C = AxB (each possibly transposed according to the flags). 
  * It is assumed that C contains enough space for the result. 
  */
 friend binary_matrix& mul(const binary_matrix& A, const bool At, const binary_matrix& B, const bool Bt, binary_matrix& C);

 void destroy() { delete[] data; memset(this,0,sizeof(binary_matrix)); }

 friend idx_t dist(const binary_matrix& A, const binary_matrix& B);

 binary_matrix& operator=(const binary_matrix& A);

private:


 /** @return the specified data block, aligned case */
 inline block_t get_block(const idx_t i, const idx_t j) const {
   return (j < last_block) ? data[i*blocks_per_row + j] : data[i*blocks_per_row + j] & trail_mask; 
 }

 /** Sets the specified data block to the given value, aligned case */
 inline void set_block(const idx_t i, const idx_t j, const block_t b)  { 
   data[i*blocks_per_row+j] = b;
 }

 /** C=A*B */
 friend binary_matrix& mul_AB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C=A^t*B */
 friend binary_matrix& mul_AtB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C=A*B^t */
 friend binary_matrix& mul_ABt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** C=A^t*B^t */
 friend binary_matrix& mul_AtBt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C);

 /** rows of the matrix */
  idx_t rows;

  /** columns of the matrix */
  idx_t cols;

  /** total number of bits. */
  idx_t len;

  /** last bit offset, may be non-zero if cols is not a multiple of block size */
  idx_t last_bit_offset; 

  /** number of actual data blocks allocated for storage */
  idx_t data_blocks;

  idx_t blocks_per_row; /** number of blocks required to hold a row */

  idx_t last_block; /** blocks_per_row-1 */

  block_t* data; 

  /** for masking the trailing bits of the matrix */
  block_t trail_mask; // mask trailing bits
};

void set_grid_width(idx_t g);

#endif
