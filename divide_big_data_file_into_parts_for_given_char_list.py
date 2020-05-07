#!/usr/bin/env python3

import string


if __name__=="__main__":

	for file in ["acc_tax.txt", "acc_tax_gb.txt"]:
	
		for ch in list(string.ascii_uppercase):

			with open(file, 'r') as in_file:

				idx = 1

				while True:
					
					line = in_file.readline().strip()
					
					idx += 1

					if line.startswith(ch):
						print("sed -n \'"+ str(idx-1)+ ", p\' "+file+" > "+ch)
						break
					
					if not line:
						break

			in_file.close()


# sed -n '2,95539889p' acc_tax.txt > A;sed -n '95539890,112443480p' acc_tax.txt > B;sed -n '112443481,218975613p' acc_tax.txt > C;sed -n '218975614,225318122p' acc_tax.txt > D;sed -n '225318123,228621271p' acc_tax.txt > E;sed -n '228621272,239634173p' acc_tax.txt > F;sed -n '239634174,297637309p' acc_tax.txt > G;sed -n '297637310,308095737p' acc_tax.txt > H;sed -n '308095738,314896924p' acc_tax.txt > I;sed -n '314896925,343266944p' acc_tax.txt > J;sed -n '343266945,351514583p' acc_tax.txt > K;sed -n '351514584,364523887p' acc_tax.txt > L;sed -n '364523888,371751005p' acc_tax.txt > M;sed -n '371751006,396957061p' acc_tax.txt > N;sed -n '396957062,441289124p' acc_tax.txt > O;sed -n '441289125,446553911p' acc_tax.txt > P;sed -n '446553912,455874967p' acc_tax.txt > Q;sed -n '455874968,459263129p' acc_tax.txt > R;sed -n '459263130,463510217p' acc_tax.txt > S;sed -n '463510218,472340232p' acc_tax.txt > U;sed -n '472340233,477204584p' acc_tax.txt > V;sed -n '477204585,$p' acc_tax.txt > W;

# sed -n '2,14651052p' acc_tax_gb.txt > Ag;sed -n '14651053,31338084p' acc_tax_gb.txt > Bg;sed -n '31338085,52393269p' acc_tax_gb.txt > Cg;sed -n '52393270,70583036p' acc_tax_gb.txt > Dg;sed -n '70583037,91478550p' acc_tax_gb.txt > Eg;sed -n '91478551,114955327p' acc_tax_gb.txt > Fg;sed -n '114955328,133763792p' acc_tax_gb.txt > Gg;sed -n '133763793,155684599p' acc_tax_gb.txt > Hg;sed -n '155684600,155755703p' acc_tax_gb.txt > Ig;sed -n '155755704,179551008p' acc_tax_gb.txt > Jg;sed -n '179551009,194218224p' acc_tax_gb.txt > Kg;sed -n '194218225,210971245p' acc_tax_gb.txt > Lg;sed -n '210971246,220421124p' acc_tax_gb.txt > Mg;sed -n '220421125,235279348p' acc_tax_gb.txt > Ng;sed -n '235279349,235379151p' acc_tax_gb.txt > Rg;sed -n '235379152,235388492p' acc_tax_gb.txt > Sg;sed -n '235388493,235485839p' acc_tax_gb.txt > Tg;sed -n '235485840,235573833p' acc_tax_gb.txt > Ug;sed -n '235573834,235575015p' acc_tax_gb.txt > Vg;sed -n '235575016,235670086p' acc_tax_gb.txt > Wg;sed -n '235670087,266484001p' acc_tax_gb.txt > Xg;sed -n '266484002,266495717p' acc_tax_gb.txt > Yg;sed -n '266495718,$p' acc_tax_gb.txt > Zg;
