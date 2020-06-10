

bool NextQuarkState(int& Spin, int& Color){
	Color++;
	if(Color>=3){
		Color=0;
		Spin++;
		if(Spin>=2){
			Spin=0;
			Color=0;
			return false;
		}
		else return true;
	}
	else return true;
}

bool AllowedCombo(int* Quark, int* Spin, int* Color, const int& NumQuarks){
	for(int iq1=0; iq1<NumQuarks; iq1++){
		for(int iq2=iq1+1; iq2<NumQuarks; iq2++){
			if(Quark[iq1]==Quark[iq2]&&Spin[iq1]==Spin[iq2]&&Color[iq1]==Color[iq2]){
				//printf("%i %i\n",iq1,iq2);
				return false;
			}
		}
	}
	return true;
}

void Quark_Blocking(){
	
	unsigned NumIter = 1000000;
	TRandom3 rangen(11);
	
	//uds anti uds
	//123 456
	
	int Proton = 112;
	int Lambda = 123;
	int Xi = 233;
	int Omega = 333;
	
	int Particle1 = Proton;
	int Particle2 = Lambda;
	
	int Quark[6];
	Quark[0] = Particle1/100;
	Quark[1] = (Particle1/10)%10;
	Quark[2] = Particle1%10;
	Quark[3] = Particle2/100;
	Quark[4] = (Particle2/10)%10;
	Quark[5] = Particle2%10;
	
	int Spin[6];
	for(int iqu=0; iqu<6; iqu++) Spin[iqu]=0;
	int Color[6];
	for(int iqu=0; iqu<6; iqu++) Color[iqu]=0;
	
	
	/*
	Spin[0] = 1;
	Spin[1] = 1;
	Spin[2] = 0;
	Spin[3] = 1;
	Spin[4] = 0;
	Spin[5] = 1;
	
	Color[0] = 2;
	Color[1] = 1;
	Color[2] = 0;
	Color[3] = 0;
	Color[4] = 1;
	Color[5] = 2;
	*/
	
	//printf("%i\n",AllowedCombo(&Quark[0],&Spin[0],&Color[0],6));
	
	//return;	
	
	
	
	
	
	//where each baryon can exist on its own
	int TotalNumStates=0;
	//where both baryons can exist simultaniously
	int AllowedNumStates=0;
	int TotalAttempts=0;
	
	for(unsigned uIter=0; uIter<NumIter; uIter++){
		for(int iqu=0; iqu<6; iqu++){
			Spin[iqu] = rangen.Integer(2);
			Color[iqu] = rangen.Integer(3);
		}
		if(AllowedCombo(&Quark[0],&Spin[0],&Color[0],3) && AllowedCombo(&Quark[3],&Spin[3],&Color[3],3))
			TotalNumStates++;
		if(AllowedCombo(&Quark[0],&Spin[0],&Color[0],6))
			AllowedNumStates++;
	}
	
	/*
	int CurrentQuark=0;
	do{
		//for(int iqu=0; iqu<6; iqu++){
	//		printf("iqu=%i\n",iqu);
	//		printf("Quark=%i\n",Quark[iqu]);
	//		printf("Spin=%i\n",Spin[iqu]);
	//		printf("Color=%i\n",Color[iqu]);
	//		printf("---------------\n");
	//	}
		//printf("\n");


		if(
		true
		//Spin[0]==0&&Spin[1]==0&&Spin[2]==0&&Spin[3]==0&&Spin[4]==0&&Spin[5]==0
		//&&
		//Color[0]==0&&Color[1]==0&&Color[2]==0&&Color[3]==0&&Color[4]==0&&Color[5]
		){
			printf("ZERO S %i\n",TotalAttempts);
			printf(" S = %i %i %i %i %i %i\n",
			Spin[0],Spin[1],Spin[2],Spin[3],Spin[4],Spin[5]);
			printf(" C = %i %i %i %i %i %i\n",
			Color[0],Color[1],Color[2],Color[3],Color[4],Color[5]);	
			usleep(500e3);
		}

		bool ChangeQuark = !(NextQuarkState(Spin[CurrentQuark],Color[CurrentQuark]));

		bool TimeToBreak = true;
		for(int iqu=0; iqu<6; iqu++){
			TimeToBreak *= Spin[iqu]==0;
			TimeToBreak *= Color[iqu]==0;
		}

		if(TimeToBreak){
			//printf("S = %i %i %i %i %i %i\n",
			//Spin[0],Spin[1],Spin[2],Spin[3],Spin[4],Spin[5]);
			//printf("C = %i %i %i %i %i %i\n",
			//Color[0],Color[1],Color[2],Color[3],Color[4],Color[5]);			
		}
		//usleep(500e3);
		//if(TimeToBreak) break;
		
		if(ChangeQuark){
			if(AllowedCombo(&Quark[0],&Spin[0],&Color[0],3) && AllowedCombo(&Quark[3],&Spin[3],&Color[3],3))
				TotalNumStates++;
			if(AllowedCombo(&Quark[0],&Spin[0],&Color[0],6))
				AllowedNumStates++;
			TotalAttempts++;
		}
		else{
			CurrentQuark++;
			if(CurrentQuark>=6){
				CurrentQuark=0;
				//break;
				//printf("HELLO\n");
			}
			else{
				
			}
		}
		//usleep(1000e3);
	}
	while(true);
	*/
	
	printf("NumIter = %u\n",NumIter);
	printf("TotalNumStates = %i\n",TotalNumStates);
	printf("AllowedNumStates = %i\n",AllowedNumStates);
	printf("Go to = %.3f\n",double(AllowedNumStates)/double(TotalNumStates));
	printf("Suppression = %.3f\n",1.-double(AllowedNumStates)/double(TotalNumStates));
	
}
