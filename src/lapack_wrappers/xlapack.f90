!> Expose lapack wrapper
module xlapack
  ! multiplication
  use vector_multiplication, only: dot_multiply, outer_product
  use general_matrix_multiplication, only: matrix_multiply
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  ! decomposition
  use singular_value_decomposition, only: svd_divide_conquer, matrix_rank
end module xlapack