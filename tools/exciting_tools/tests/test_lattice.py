"""
Tests for functions in lattice.py 
"""
import numpy as np

from excitingtools.lattice import reciprocal_lattice_vectors, parallelpiped_volume, plane_transformation
from excitingtools.math_utils import triple_product


def test_parallelpiped_volume():

    # FCC lattice vectors with arbitrary lattice constant
    a = 3.2 * np.array([[0., 1., 1.],
                        [1., 0., 1.],
                        [1., 1., 0.]])

    triple_product_result = triple_product(a[:, 0], a[:, 1], a[:, 2])
    volume = parallelpiped_volume(np.array([a[:, 0], a[:, 1], a[:, 2]]))

    assert np.isclose(triple_product_result, 65.5360)
    assert np.isclose(volume, 65.5360)

    triple_product_result = triple_product(a[:, 0], a[:, 2], a[:, 1])
    volume = parallelpiped_volume(np.array([a[:, 0], a[:, 2], a[:, 1]]))

    assert np.isclose(triple_product_result, -65.5360), "triple product can be negative"
    assert np.isclose(volume, 65.5360), "volume is always defined as positve"


def test_reciprocal_lattice_vectors():
    r"""
    Test
        \mathbf{a}_i \cdot \mathbf{b}_j = 2 \pi \delta_{ij}

    """

    # fcc with arbitrary lattice constant
    a = 3.2 * np.array([[0., 1., 1.],
                        [1., 0., 1.],
                        [1., 1., 0.]])
    b = reciprocal_lattice_vectors(a)

    assert a.shape == (3, 3)
    assert b.shape == (3, 3)

    a_dot_b = np.transpose(a) @ b
    off_diagonal_indices = np.where(~np.eye(a_dot_b.shape[0], dtype=bool))
    off_diagonal_elements = a_dot_b[off_diagonal_indices]

    assert np.allclose(a_dot_b.diagonal(), 2 * np.pi)
    assert np.allclose(off_diagonal_elements, 0.)

def test_plane_transformation():
    plot_vec = np.array([[-0.042506, -0.042506,  0.050235],
                         [-0.042506,  0.050235, -0.042506],
                         [-0.085013,  0.007728,  0.007728]]
                       )

    rec_lat_vec = np.array([[-0.5881478,  0.5881478,  0.5881478],
                            [ 0.5881478, -0.5881478,  0.5881478],
                            [ 0.5881478,  0.5881478, -0.5881478]]
                          )

    transformation_matrix = plane_transformation(plot_vec, rec_lat_vec) # This is supposed to be 3x3 Identity

    diagonal = transformation_matrix.diagonal()
    offdiagonal = transformation_matrix[np.where(~np.eye(transformation_matrix.shape[0], dtype=bool))]

    vec=np.array([-4.545696326767e-3, 4.999971885544e-2, 5.881478006738e-7])

    assert np.allclose(transformation_matrix.dot(plot_vec)[2,:], 0.0, atol=1e-6)

def test_plane_transformation():
    plot_vec = np.array([[-0.042506, -0.042506,  0.050235],
                         [-0.042506,  0.050235, -0.042506],
                         [-0.085013,  0.007728,  0.007728]]
                       )

    rec_lat_vec = np.array([[-0.5881478,  0.5881478,  0.5881478],
                            [ 0.5881478, -0.5881478,  0.5881478],
                            [ 0.5881478,  0.5881478, -0.5881478]]
                          )

    transformation_matrix = plane_transformation(plot_vec, rec_lat_vec) # This is supposed to be 3x3 Identity

    diagonal = transformation_matrix.diagonal()
    offdiagonal = transformation_matrix[np.where(~np.eye(transformation_matrix.shape[0], dtype=bool))]

    vec=np.array([-4.545696326767e-3, 4.999971885544e-2, 5.881478006738e-7])

    assert np.allclose(transformation_matrix.dot(plot_vec)[2,:], 0.0, atol=1e-6)