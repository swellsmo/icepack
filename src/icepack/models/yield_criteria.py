


def vonMises(self, **kwargs):
  M = itemgetter("membrane_stress")(**kwargs)
  σ_e = sqrt(inner(M, M) - det(M))
  return σ_e 
