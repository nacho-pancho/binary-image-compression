#include "binmat.h"
#include <iomanip>
#include <cstring>
#include <bitset>
#include <cassert>
#include <omp.h>


/* information on how to vectorize these operations using x86 SSE extensions in C is
   described in https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html */

/* https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable */

/*
 * very efficient general purpose method for counting bits w/o requiring
 * a look up table. An alternative is a LUT as a 32 bit block where one has 16 2 bit
 * counts for every of the 16 possible 4 bit nibbles.
 * 
 * There is a POPCNT CPU extension in intel-base 64 bit processors since 2014,
 * but I don't want to be compiler dependent for now.
 */
static idx_t block_weight(const block_t& v)
{
  static const unsigned char byte_weight_table[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};
  const unsigned char * p = (const unsigned char *) &v;
  register idx_t w = 0;
  for (idx_t i = 0; i < sizeof(block_t); ++i) {
    w += byte_weight_table[p[i]];
  }
  return w;
}

/* parity computation, much faster than sum */
static bool block_sum(const block_t& v) {
  static const bool byte_parity_table[256] = {
#   define P2(n) n, n^1, n^1, n
#   define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#   define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
    P6(0), P6(1), P6(1), P6(0)
  };
  const unsigned char * p = (unsigned char *) &v;
  idx_t ti = 0;
#pragma omp parallel for
  for (idx_t i = 0; i < sizeof(block_t); ++i) {
    ti ^= p[i];
  }
  return byte_parity_table[ti];
  //  return ParityTable256[tip[0] ^ p[1] ^ p[2] ^ p[3]];
}

idx_t binary_matrix::weight() const {
  if (rows*cols == 0)
    return 0;
  idx_t w = 0;
  for (idx_t i = 0; i < rows; ++i) {
    for (idx_t j = 0; j < blocks_per_row ; ++j) {
      w += block_weight(get_block(i,j));
    }
  }
  return w;
}

idx_t binary_matrix::row_weight(idx_t i) const {
  //  assert(i >= 0);
  assert(i < rows);
  if (rows == 0) return 0;
  idx_t w = 0;
  for (idx_t j = 0; j < blocks_per_row; ++j) {
    w += block_weight(get_block(i,j));
  }
  return w;
}

idx_t binary_matrix::col_weight(idx_t j) const {
  //  assert(j >= 0);
  assert(j < cols);
  if (cols == 0) return 0;
  block_t* pd = &data[j / BITS_PER_BLOCK];
  const block_t mask = MSB >> (j % BITS_PER_BLOCK);
  idx_t w = 0;
  for (idx_t i = 0; i < rows; i+=blocks_per_row) {
    if (pd[i] & mask)
      w++;
  }
  return w;
}

bool binary_matrix::sum() const {
  if (rows*cols == 0)
    return 0;
  idx_t w = 0;
  for (idx_t i = 0; i < rows; ++i) {
    for (idx_t j = 0; j < blocks_per_row; ++j) {
      w ^= block_sum(get_block(i,j));
    }
  }
  return w;
}


bool binary_matrix::row_sum(idx_t i) const {
  //  assert(i >= 0);
  assert(i < rows);
  if (cols == 0) return 0;
  idx_t w = 0;
  for (idx_t j = 0; j < blocks_per_row; ++j) {
    w ^= block_sum(get_block(i,j));
  } 
  return w;
}

bool binary_matrix::col_sum(idx_t j) const {
  block_t* pd = &data[j / BITS_PER_BLOCK];
  block_t mask = MSB >> (j % BITS_PER_BLOCK);
  bool w = false;
  for (idx_t i = 0; i < data_blocks; i+=blocks_per_row) {
    if (pd[i] & mask) { w = !w; }
  }
  return w;
}

binary_matrix::binary_matrix(const binary_matrix& other)
  :rows(other.rows), cols(other.cols), len(other.len)
{ 
  blocks_per_row = other.blocks_per_row;
  data_blocks = other.data_blocks;
  last_bit_offset = other.last_bit_offset;
  trail_mask = other.trail_mask;
  last_block = other.last_block;
  data = new block_t[data_blocks];
  std::copy(other.data,other.data+other.data_blocks,data);
}

binary_matrix::binary_matrix(idx_t _rows, idx_t _cols)
  :rows(_rows), cols(_cols), len(_rows*_cols)
{ 
  blocks_per_row = (cols+BITS_PER_BLOCK-1)/BITS_PER_BLOCK; // ceil
  data_blocks = blocks_per_row*rows;
  data = new block_t[data_blocks];
  last_bit_offset = (cols-1) % BITS_PER_BLOCK;
  trail_mask = ONES << (BITS_PER_BLOCK-last_bit_offset-1);
  last_block = blocks_per_row - 1;
}

void binary_matrix::allocate(idx_t _rows, idx_t _cols) {
  rows = _rows;
  cols = _cols;
  len = rows*cols;
  blocks_per_row = (cols+BITS_PER_BLOCK-1)/BITS_PER_BLOCK; // ceil
  data_blocks = blocks_per_row*rows;
  data = new block_t[data_blocks];
  last_bit_offset = (cols-1) % BITS_PER_BLOCK;
  trail_mask = ONES << (BITS_PER_BLOCK-last_bit_offset-1);
  last_block = blocks_per_row - 1;
}

void binary_matrix::clear() { 
  if (rows*cols == 0) return;
  memset(data,0,sizeof(block_t)*rows*blocks_per_row); 
}

void binary_matrix::set() {
  if (rows*cols == 0) return;
  memset(data,0xff,sizeof(block_t)*rows*blocks_per_row); 
}

void binary_matrix::flip() {
  if (rows*cols == 0) return;
  for (idx_t i = 0; i < data_blocks; ++i) {
    data[i] = ~data[i];
  }
}

binary_matrix& binary_matrix::operator=(const binary_matrix& A) {
   delete[] data;
   memcpy(this,&A,sizeof(binary_matrix));
   return *this;
}


binary_matrix binary_matrix::get_copy() const {
  binary_matrix A = *this; // shallow copy
  // now replace memory
  A.data = new block_t[A.data_blocks];
  memcpy(A.data,data,sizeof(block_t)*data_blocks);
  return A;
}

void binary_matrix::copy_to(binary_matrix& A) const {
  memcpy(A.data,data,sizeof(block_t)*data_blocks);
}

void binary_matrix::transpose_to(binary_matrix& A) const {
  assert(rows == A.cols);
  assert(cols == A.rows);
  const idx_t n = cols;
  binary_matrix v(1,rows);
  for (idx_t j = 0; j < n; ++j) {
    copy_col_to(j,v);
    A.set_row(j,v);
  }
}

binary_matrix binary_matrix::get_transposed() const {
  binary_matrix A(rows,cols);
  transpose_to(A);
  return A;  
}

binary_matrix binary_matrix::get_col(const idx_t j) const {
  //  assert(0 <= j);
  assert(j <= cols);
  binary_matrix col(1,rows); // NOTE: it is a ROW vector
  copy_col_to(j,col);
  return col;
}

void binary_matrix::copy_col_to(const idx_t j, binary_matrix& col) const {
  //  assert(0 <= j);
  assert(j <= cols);
  col.clear();
  const block_t mask1 = MSB >> (j % BITS_PER_BLOCK);
  block_t mask2 = MSB;
  block_t* const pd = data + j/BITS_PER_BLOCK;
  for (idx_t i = 0, i2 = 0; i < rows; ++i) {
    if (pd[i*blocks_per_row] & mask1) 
      col.data[i2] |= mask2;
    mask2 >>=1;
    if (!mask2) {
      mask2 = MSB;
      ++i2;
    }
  }
}

binary_matrix binary_matrix::get_row(const idx_t i) const {
  //  assert(0 <= i);
  assert(i < rows);
  binary_matrix row(1,cols);
  copy_row_to(i,row);
  return row;
}

void binary_matrix::copy_row_to(const idx_t i, binary_matrix& row) const {
  //  assert(0 <= i);
  assert(i < rows);
  block_t* const pd = data+i*blocks_per_row;  
  for (idx_t j = 0; j < blocks_per_row; ++j) {
    row.data[j] = pd[j];
  }
}

binary_matrix binary_matrix::get_submatrix(const idx_t i0, const idx_t i1, const idx_t j0, const idx_t j1) const {
  assert(i0 < i1);
  assert(j0 < j1);
  binary_matrix B(i1-i0,j1-j0);
  copy_submatrix_to(i0,i1,j0,j1,B);
  return B;
}

void  binary_matrix::copy_submatrix_to(const idx_t i0, const idx_t i1, const idx_t j0, const idx_t j1, binary_matrix& B) const {
  // could remove these in the future
  assert(i0 < i1);
  assert(j0 < j1);
  //B.set(); // DEBUG
  const idx_t boff = j0 % BITS_PER_BLOCK;
  const idx_t doff = j0 / BITS_PER_BLOCK + i0*blocks_per_row;
  //  std::cout << "get_sub: i0=" << i0 << " j0=" << j0 << " boff=" << boff << " doff=" << doff  << std::endl;
  if (boff == 0) { // submatrix is aligned to word boundary
    for (idx_t di = 0; di < B.rows; di++) {
      for (idx_t dj = 0; dj < B.blocks_per_row; dj++) {
	idx_t k = doff + di*this->blocks_per_row + dj;	
	idx_t sb = (k < data_blocks)? data[k] : 0;
	B.set_block(di,dj,sb);
      }
    }
  } else { // not aligned: a little more difficult
    for (idx_t di = 0; di < B.rows; di++) {
      for (idx_t dj = 0; dj < B.blocks_per_row; dj++) {
	const idx_t k1 = doff + di*this->blocks_per_row + dj;
	const idx_t k2 = k1 + 1;
	const idx_t s1 = (k1 < data_blocks)? data[k1] : 0;
	//	if (k2 < data_blocks) {	
	  const idx_t s2 = (k2 < data_blocks)? data[k2] : 0;
	  B.set_block( di,dj,(s1 << boff) | (s2 >> (BITS_PER_BLOCK-boff)));
	  //	} else {
	  //	  B.set_block( di,dj,(s1 << boff));
	  //}
      }
    }
  }
}

binary_matrix binary_matrix::get_vectorized() const {
  binary_matrix v(1,rows*cols);
  copy_vectorized_to(v);
  return v;
}

void binary_matrix::copy_vectorized_to(binary_matrix& v) const {
  v.clear();
  for (idx_t si = 0; si < rows; ++si) {
    const idx_t loff = this->cols*si;
    const idx_t boff = loff % BITS_PER_BLOCK;
    idx_t dj = loff / BITS_PER_BLOCK;
    for (idx_t sj = 0; sj < blocks_per_row; sj++,dj++) {
      block_t sb = this->get_block(si,sj);
      v.data[dj] |= (sb >> boff);
      if (dj < (v.data_blocks-1)) {
	v.data[dj+1] = (sb << (BITS_PER_BLOCK-boff));
      }
    }
  }
}

void binary_matrix::set_vectorized(const binary_matrix& src)
{
  clear();
  for (idx_t di = 0; di < rows; ++di) {
    const idx_t loff = this->cols*di; // same for each dest row
    const idx_t boff = loff % BITS_PER_BLOCK; // same for each destination row
    idx_t sj = loff / BITS_PER_BLOCK;
    //    std::cout << "loff=" << loff << " boff=" << boff << " sj=" << sj << std::endl;
    //    std::cout << *this << std::endl;
    // BUG HERE: reading ony byte past the end of src!
    // 
    for (idx_t dj = 0; dj < blocks_per_row; dj++,sj++) {
      idx_t k = di*blocks_per_row + dj;
      data[k]  |= (src.data[sj] << boff);
      if (boff && (sj < (src.data_blocks-1))) {  
	data[k]  |= (src.data[sj+1] >> (BITS_PER_BLOCK-boff));
      }
    }
  }
}

void binary_matrix::set_col(const idx_t j, const binary_matrix& B)
{
  //  assert(0 <= j);
  assert(j <= cols);
  block_t masks = MSB;
  const block_t maskd = MSB >> (j % BITS_PER_BLOCK);
  block_t* const pd = data + j/BITS_PER_BLOCK;
  for (idx_t id = 0, is = 0; id < rows; ++id) {
    pd[id*blocks_per_row] &= ~maskd;
    if (B.data[is] & masks) 
      pd[id*blocks_per_row] |= maskd;
    masks >>=1;
    if (!masks) {
      masks = MSB;
      ++is;
    }
  }
}

void binary_matrix::set_row(const idx_t i, const binary_matrix& B) 
{
  //  assert(0 <= i);
  assert(i < rows);
  block_t* const pd = data+i*blocks_per_row;  
  block_t* const ps = B.data;  
  for (idx_t j = 0; j < blocks_per_row; ++j) {
    pd[j] = ps[j];
  }
}

void binary_matrix::set_submatrix(const idx_t i0,const  idx_t j0, const binary_matrix& B) 
{
  // could remove these in the future
  const idx_t boff = j0 % BITS_PER_BLOCK;
  const idx_t lastoff = ( (boff+B.cols) % BITS_PER_BLOCK );
  const idx_t spanned_blocks = (BITS_PER_BLOCK-1 + boff + B.cols) / BITS_PER_BLOCK;
  if ((boff == 0) || (spanned_blocks==1)) { // easy cases
    const block_t maskr = boff ? ~(ONES<<(BITS_PER_BLOCK-boff)) : ONES;
    const block_t maskl = lastoff ? (ONES <<(BITS_PER_BLOCK - lastoff)) : ONES;
    const block_t mask = maskr & maskl;
    //    std::cout << "set_sub (a): i0=" << i0 << " j0=" << j0 << " boff=" << boff << " lastoff=" << lastoff << " mask=" << bm_bitset(mask) << std::endl;
    for (idx_t si = 0,di = i0; (si < B.rows) && (di < this->rows); si++,di++) {
      idx_t dj = j0/BITS_PER_BLOCK;
      idx_t sj = 0;
      for ( ; (sj < B.last_block) && (dj < this->blocks_per_row); sj++,dj++) {
	this->set_block(di,dj, B.get_block(si,sj));
      }
      this->set_block(di,dj, (get_block(di,dj) & ~mask) | ((B.get_block(si,sj)>>boff) & mask));
    }
  } else { // more difficult! Each block in B covers more than one block in this matrix
    const block_t mask2(ONES << (BITS_PER_BLOCK-boff));
    const block_t mask1 = ~mask2;
    const block_t mask3(ONES << (BITS_PER_BLOCK-lastoff));
    //    std::cout << "set_sub (na): i0=" << i0 << " j0=" << j0 << " boff=" << boff << " spanned=" << spanned_blocks << " mask1=" << bm_bitset(mask1) << " mask2="<< bm_bitset(mask2) << " mask3="<< bm_bitset(mask3) << std::endl;
    for (idx_t si = 0, di=i0; (si < B.rows) && (di < this->rows); si++,di++) {
      idx_t dj = j0 / BITS_PER_BLOCK;
      idx_t sj = 0;
      for ( ; (sj < B.last_block) && (dj < this->last_block); sj++,dj++) {
        block_t sb = B.get_block(si,sj);
	this->set_block(di,dj,  (this->get_block(di,dj) & ~mask1) | (sb >> boff));
        this->set_block(di,dj+1, (this->get_block(di,dj+1) & ~mask2) | (sb<<(BITS_PER_BLOCK-boff)));
      }
      // last block: second part of block may be shorter than boff
      block_t sb = B.get_block(si,sj);
      block_t db = (dj <= last_block) ? get_block(di,dj) : 0;
      set_block(di,dj, (db & ~mask1) | (sb >> boff));
      if (dj < this->last_block) {
        this->set_block( di,dj+1, (this->get_block(di,dj+1) & ~mask3) | ((sb << (BITS_PER_BLOCK-boff)) & mask3) ); 
      }
    }
  }
}


void binary_matrix::add_rows(idx_t nrows)  {
  rows += nrows;
  idx_t new_data_blocks = blocks_per_row*rows;
  block_t* new_data = new block_t[new_data_blocks];
  std::copy(data,data+data_blocks,new_data);
  data_blocks = new_data_blocks;
  delete[] data;
  data = new_data;
}
void binary_matrix::remove_rows(idx_t nrows)  {
  if (nrows < rows) {
    rows -= nrows;  
    data_blocks -= blocks_per_row*nrows;
  } else {
    data_blocks = 0;
    rows = 0;
  }
}
#if 0

void binary_matrix::add_row(binary_matrix& v)  {
  rows ++;
  idx_t new_data_blocks += blocks_per_row;
  block_t* new_data = new block_t[new_data_blocks];
  std::copy(data,data+data_blocks,new_data);
  data_blocks = new_data_blocks;
  delete[] data;
  data = new_data;
  set_row(rows-1,v);
}

void binary_matrix::remove_row(idx_t i)  {
  assert(i < rows);
  assert(i >= 0);
  if (nrows < rows} {
    std::copy(data+blocks_per_row*(i+1),data+data_blocks,data+blocks_per_row*i);
    rows -= nrows;  
    data_blocks -= blocks_per_row*nrows;
  } else {
    data_blocks = 0;
    rows = 0;
  }
}
#endif

// C is assumed to have been allocated and have the appropriate dimension
binary_matrix& add(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.rows == B.rows);
  assert(C.cols == A.cols);
  assert(C.cols == B.cols);
  const idx_t M = A.rows;
  const idx_t Nb = C.blocks_per_row;
  for (idx_t i = 0; i < M; ++i) {      
    for (idx_t j = 0; j < Nb; ++j) {
      C.set_block(i,j, A.get_block(i,j) ^ B.get_block(i,j));
    }
  } 
  return C;
}

// C is assumed to have been allocated and have the appropriate dimension
binary_matrix& bool_and(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.rows == B.rows);
  assert(C.cols == A.cols);
  assert(C.cols == B.cols);
  const idx_t M = A.rows;
  const idx_t Nb = C.blocks_per_row;
  for (idx_t i = 0; i < M; ++i) {      
    for (idx_t j = 0; j < Nb; ++j) {
      C.set_block(i,j, A.get_block(i,j) & B.get_block(i,j));
    }
  } 
  return C;
}


idx_t dist(const binary_matrix& A, const binary_matrix& B)
{
  assert(A.rows == B.rows);
  assert(A.cols == B.cols);
  const idx_t M = A.rows;
  const idx_t Nb = A.blocks_per_row;
  idx_t w = 0;
  for (idx_t i = 0; i < M; ++i) {      
    for (idx_t j = 0; j < Nb; ++j) {
      w += block_weight(A.get_block(i,j) ^ B.get_block(i,j));
    }
  } 
  return w;
}


// C is assumed to have been allocated and have the appropriate dimension
binary_matrix& mul_AB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.cols == B.cols);
  assert(A.cols == B.rows);
  const idx_t M = A.rows;
  const idx_t Ns = B.blocks_per_row;
  const idx_t K = B.rows;
  C.clear();
  block_t mask = MSB;
  block_t koff = 0;
  for (idx_t k = 0; k < K; ++k) {
    for (idx_t i = 0; i < M; ++i) {
      if (A.get_block(i,koff) & mask) { // if nonzero, the i-th row of C is updated by adding the k-th row of B
	for (idx_t j = 0; j < Ns; ++j) {
	  C.set_block(i,j, C.get_block(i,j) ^ B.get_block(k,j));
	}
      } 
    }
    mask >>= 1;
    if (mask == 0) {
      mask = MSB;
      koff++; 
    }
  }
  return C;
}

binary_matrix& mul_AtB(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.cols);
  assert(C.cols == B.cols);
  assert(A.rows == B.rows);
  const idx_t K = A.rows;
  const idx_t M = A.cols;
  const idx_t Ns = B.blocks_per_row;
  C.clear();
  for (idx_t k = 0; k < K; ++k) { 
    block_t mask = MSB;
    block_t ioff = 0;
    for (idx_t i = 0; i < M; ++i) {      
      //      std::cout << "i=" << i << "mask=" << std::bitset<BITS_PER_BLOCK>(mask) << " C=" << C << std::endl;   
      if (A.get_block(k,ioff) & mask) { // if nonzero, the i-th row of C is updated by adding the k-th row of B
	for (idx_t j = 0; j < Ns; ++j) {
	  C.set_block(i,j, C.get_block(i,j) ^ B.get_block(k,j));
	}
      } 
      mask >>= 1;
      if (mask == 0) {
	mask = MSB;
	ioff++;
      }
    } // i
  } // k
  return C;
}

binary_matrix& mul_ABt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.rows);
  assert(C.cols == B.rows);
  assert(A.cols == B.cols);
  const idx_t M = A.rows;
  const idx_t N = B.cols;
  const idx_t K = A.blocks_per_row;
  for (idx_t i = 0; i < M; ++i) {
    for (idx_t j = 0; j < N; ++j) {
      bool a = 0;
      for (idx_t k = 0; k < K; ++k) {
	a ^= block_sum(A.get_block(i,k) & B.get_block(j,k)); // pA[k] & pB[k]);
      }
      C.set(i,j,a);
    }
  }
  return C;
}

binary_matrix& mul_AtBt(const binary_matrix& A, const binary_matrix& B, binary_matrix& C)
{
  assert(C.data != 0);
  assert(C.rows == A.cols);
  assert(C.cols == B.rows);
  assert(A.rows == B.cols);
  // FALTA!
  return C;
}

binary_matrix& mul(const binary_matrix& A, const bool At, const binary_matrix& B, const bool Bt, binary_matrix& C) {
  if (At && Bt) {
    return mul_AtBt(A,B,C);
  } else if (At && !Bt) {
    return mul_AtB(A,B,C);
  } else if (!At && Bt) {
    return mul_ABt(A,B,C);
  } else {
    return mul_AB(A,B,C);
  }
}

static idx_t grid_width = 10;

void set_grid_width(idx_t g) { grid_width = g; }

// slow: I don't think anyone cares about fast dumping, since most of the time
// will be I/O anyway.
std::ostream& operator<<(std::ostream& out, const binary_matrix& A)  {
  out << "rows=" << A.rows << "\tcols=" << A.cols << "\tlen=" << A.len 
      << "\tbpw="<< BITS_PER_BLOCK 
      << "\tdw=" << A.data_blocks
      << "\twpr=" << A.blocks_per_row
      << "\ttm=" << bm_bitset(A.trail_mask);
  out << std::endl;
  out << "       ";
  for (idx_t j = 0; j < A.cols; ++j)
    out << ((j % BITS_PER_BLOCK) ? ' ' : '|') << ' ';
  out << std::endl;
  for (idx_t i = 0; i < A.rows; ++i) {
    out << std::setw(5)<< i << "  "; // << ((i % grid_width) ? " | " : " + ");
    for (idx_t j = 0; j < A.cols; ++j) {
      char pixel = (i % grid_width) || (j % grid_width) ? '.' : '+';
      out << std::setw(1) << (A.get(i,j) ? '#':pixel) << ' ';
    }
    out << std::endl;
  }
  return out;
}

