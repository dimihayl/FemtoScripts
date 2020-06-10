{
	const double Mass1 = 938.272;
	//const double Mass1 = 139.57;
	const double Mass2 = 1321.71;
	//const double Mass2a = 1115.683;
	const double RelMom = 110;//k*
	
	//const double MASS = sqrt(Mass1*Mass1+Mass2*Mass2+2*sqrt(Mass1*Mass1+RelMom*RelMom)*sqrt(Mass2*Mass2+RelMom*RelMom)+2*RelMom*RelMom);
	const double MASS = sqrt(Mass1*Mass1+RelMom*RelMom)+sqrt(Mass2*Mass2+RelMom*RelMom);
	
	printf("RelMom = %.4f\n",RelMom);
	printf("MASS = %.4f\n",MASS);
	printf("THR = %.4f\n",Mass1+Mass2);
	printf("E = %.4f\n",MASS-Mass1-Mass2);
	//printf("Ea = %.4f\n",MASS-Mass2a-Mass2a);
}
