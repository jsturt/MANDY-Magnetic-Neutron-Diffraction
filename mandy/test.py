import mandyCrystal as mC


crys = mC.mandyCrystal(r'C:\Users\John\Documents\Uni\2022 Internship\Python Code\MANDY-Magnetic-Neutron-Diffraction-off-axis-moments\NbFe2 Example\NbFe2_mp-568901_primitive.cif')
crys.build()
crys.createTransverseSDW(0.1)

#values = [row3[1:4] for row3 in crys.pos_df.itertuples() ]


