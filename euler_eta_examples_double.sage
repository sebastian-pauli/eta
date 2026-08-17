# Copyright 2026 Sebastian Pauli
#
# This program is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 
# of the License, or any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program. If not, see
# <https://www.gnu.org/licenses/>.

attach("euler_eta.sage")

print("Find zero s0 of eta(s)-1 connected to the zero s=-6 of eta")
L0 = euler_eta_c_path(mpc(-6))
s0 = L0[-1][1]	
print("s0 =",s0)

print("Find first approximation s1 of double zero of a fractional derivative")
L1 = euler_eta_frac_path(s0,max_real=-4,delta=10^-2)
for bs in L1:
  if bs[1].imag < 10^-2: break
a1 = bs[0]; s1 = bs[1]
print("approximate order of derivative a1 =",a1)
print("near s1 =",s1)
print("a1-th derivative of euler_eta(s1) is",euler_eta(s1,a1))
print("(a1+1)-th derivative of euler_eta(s1) is",euler_eta(s1,a1+1))


print("Improve approximation of double zero")
sb = findroot([lambda s,alpha: euler_eta(s,alpha),lambda s,alpha: euler_eta(s,alpha+1)],(s1.real,mpf(a1)))
s2 = sb.tolist()[0][0]; a2 = sb.tolist()[1][0]
alpha = a2.real; z = s2.real
print("Double zeros is near")
print("approximate order of derivative alpha =",alpha)
print("near z =",z)
print("alpha-th derivative of euler_eta(s3) is",euler_eta(z,alpha))
print("(alpha+1)-th derivative of euler_eta(s3) is",euler_eta(z,alpha))

G = point([(RR(s[1].real),RR(s[1].imag)) for s in L0], color="gray",legend_label="$\\eta(s)-c=0$ for $c\\in[0,1]$")
G+= point([(RR(s[1].real),RR(s[1].imag)) for s in L1], color="cyan",legend_label="$\\eta^{(\\alpha)}(\\sigma+it)=0$")
G+= point((RR(z),0),color="red",legend_label="double zero")
G.show()

print("Plot of alpha-th and (alpha+1)-th derivative of eta on [-4.9,-4.85]")
p0r = plot(lambda x:RR(euler_eta(x,alpha).real),-4.9,-4.85,color="blue",legend_label="$\\Re(\\eta^{(\\alpha)})(\\sigma)$")
p0i = plot(lambda x:RR(euler_eta(x,alpha).imag),-4.9,-4.85,color="red",legend_label="$\\Im(\eta^{(\\alpha)})(\\sigma)$")
p1r = plot(lambda x:RR(euler_eta(x,alpha+1).real),-4.9,-4.85,color="cyan",legend_label="$\\Re(\\eta^{(\\alpha+1)})(\\sigma)$")
p1i = plot(lambda x:RR(euler_eta(x,alpha+1).imag),-4.9,-4.85,color="magenta",legend_label="$\Im(\\eta^{(\\alpha+1)})(\\sigma)$")
(p0r+p0i+p1r+p1i).show()

print("Count number of zeros with argument principle")
d = 10**-6
nr0 = euler_eta_number_of_zeros(z-d-d*I,z+d+d*I,alpha)
print("The number of zeros (counting multiplicity) in the box with sidelength 2*10^-6 around z is",nr0)
nr1 = euler_eta_number_of_zeros(z-d-d*I,z+d+d*I,alpha+1)
print("The number of zeros (counting multiplicity) in the box with sidelength 2*10^-6 around z is",nr1)

