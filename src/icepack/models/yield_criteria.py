# Copyright (C) 2026 by Sarah Wells-Moran <swellsmo@uchicago.edu>
#
# This file is modified from icepack.
#
# icepack is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The full text of the license can be found in the file LICENSE in the
# icepack source directory or at <http://www.gnu.org/licenses/>.

from operator import itemgetter
import firedrake
from firedrake import inner, sqrt, det, min_value, max_value, conditional
from .transport import TransportEquation
from icepack.constants import year
from icepack.utilities import eigenvalues


def vonMises(M):
  r"""Compute the von Mises stress for a given membrane stress"""
  σ_e = sqrt(inner(M, M) - det(M))
  return σ_e 

def tresca(M):
  r"""Compute the Tresca/maximum shear stress for a given membrane stress"""
  σ1, σ2 = eigenvalues(M)
  σ_e = max_value(abs(σ1), abs(σ2), abs(σ1-σ2))
  return σ_e 

def mohrCoulomb(M):
  r"""Compute the von Mises stress for a given membrane stress"""
  σ1, σ2 = eigenvalues(M)
  σ_e = sqrt(inner(M, M) - det(M))
  return σ_e 

def druckerPrager(M):
  r"""Compute the von Mises stress for a given membrane stress"""
  σ1, σ2 = eigenvalues(M)
  σ_e = sqrt(inner(M, M) - det(M))
  return σ_e 

def hayhurst(M, **kwargs):
  r"""Compute the von Mises stress for a given membrane stress"""
  alpha = kwargs.get("alpha", 0.21)
  beta = kwargs.get("beta", 0.63)
  σ1, σ2 = eigenvalues(M)
  σ_e = sqrt(inner(M, M) - det(M))
  return σ_e 
