{
	const double Mass1 = 938.272;
//p = 938.272
//n = 939.565
//Lambda = 1115.683
//sigma 0 = 1192.642
//xi- = 1321.71
	//const double Mass1 = 139.57;
	const double Mass2 = 1115.683;
	//const double Mass2a = 1115.683;
	//const double ERES = 37;
	//const double MASS = Mass2a+Mass2a+ERES;
	const double MASS = 939.565+1192.642;	

	const double RelMom = sqrt(MASS*MASS*0.25+pow((Mass1*Mass1-Mass2*Mass2)*0.5/MASS,2.)-(Mass1*Mass1+Mass2*Mass2)*0.5);
	
	printf("RelMom = %.4f\n",RelMom);
	printf("MASS = %.4f\n",MASS);
	printf("THR = %.4f\n",Mass1+Mass2);
	printf("E = %.4f\n",MASS-Mass1-Mass2);
	//printf("Ea = %.4f\n",MASS-Mass2a-Mass2a);
}
