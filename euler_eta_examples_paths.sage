# Copyright 2026 by Sebastian Pauli
#
# This program is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 
# of the License, or any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <https://www.gnu.org/licenses/>.

attach("euler_eta.sage")

# zeros of Riemann zeta function are used as starting points
zetazeros = [ 14.1347251417346937904572519835625, 21.0220396387715549926284795938969, 25.0108575801456887632137909925628, 30.4248761258595132103118975305840, 32.9350615877391896906623689640747, 37.5861781588256712572177634807053, 40.9187190121474951873981269146334, 43.3270732809149995194961221654068, 48.0051508811671597279424727494277, 49.7738324776723021819167846785638,
 49.7738324776723021819167846785638,52.9703214777144606441472966088808,56.4462476970633948043677594767060,59.3470440026023530796536486749922,60.8317785246098098442599018245241]

print("Plot of path starting at first non trivial zero of zeta")
print("eta-c path")
s = mpc(1/2+I*zetazeros[0])
print("start at",s)
L0 = euler_eta_c_path(s)
print("fractional derivative path")
s = L0[-1][1]
print("start at",s)
L1 = euler_eta_frac_path(s)
G = point([(RR(s[1].real),RR(s[1].imag)) for s in L0], color="gray",legend_label="$\\eta(s)-c=0$ for $c\\in[0,1]$")
G+= point([(RR(s[1].real),RR(s[1].imag)) for s in L1], color="cyan",legend_label="$\\eta^{(\\alpha)}(\\sigma+it)=0$")
G+= plot_labels(L1)
G.show()

print("Adding plots of more paths, this takes a while")
for z in zetazeros[1:]:
  s = mpc(1/2+I*z)
  print("eta-c path")
  print("start at",s)
  L0 = euler_eta_c_path(s)
  print("fractional derivative path")
  s = L0[-1][1]
  print("start at",s)
  L1 = euler_eta_frac_path(s)
  G+= point([(RR(s[1].real),RR(s[1].imag)) for s in L0], color="gray")
  G+= point([(RR(s[1].real),RR(s[1].imag)) for s in L1], color="cyan")
  G+= plot_labels(L1)
G.show(figsize=(8,12))

print("Adding plots of paths starting at zeros z with Re(z)=1")
for n in range(1,7):
  s = mpc(1+2*n*pi()*I/log(2))
  print("start at",s)
  L0 = euler_eta_c_path(s)
  G+=point([(RR(s[1].real),RR(s[1].imag)) for s in L0], color="gray")
G.show(figsize=(8,12))
