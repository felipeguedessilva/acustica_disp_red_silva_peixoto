import numpy as np

# Define fuzzy sets
a = np.array([0.2, 0.5, 0.7, 1.0])
b = np.array([0.1, 0.3, 0.9, 0.6])

# Fuzzy union
c = np.fmax(a, b)
print("Fuzzy Union:", c)

# Fuzzy intersection
d = np.fmin(a, b)
print("Fuzzy Intersection:", d)

# Fuzzy complement
e = 1 - a
print("Fuzzy Complement:", e)
