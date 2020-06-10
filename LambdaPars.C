void LambdaPars(){

    const bool Compact = true;

    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feed-down (flat)
    //3 = misidentified
    const unsigned NumChannels_p = 4;
    const double ProtonPurity = 0.99;
    const double Purities_p[NumChannels_p] = {ProtonPurity,ProtonPurity,ProtonPurity,1.-ProtonPurity};
    const double Fraction_p[NumChannels_p] = {0.87,0.13*0.7,0.13*0.3,1};

    //for the Lambda:
    //0 = primary
    //1 = from Sigma0
    //2 = from Xim
    //3 = from Xi0
    //4 = misidentified
    const unsigned NumChannels_L = 5;
    const double LambdaPurity = 0.96;
    const double Purities_L[NumChannels_L] = {LambdaPurity,LambdaPurity,LambdaPurity,LambdaPurity,1.-LambdaPurity};
    const double Fraction_L[NumChannels_L] = {0.76*3./4.,0.76*1./4.,0.24*0.5,0.24*0.5,1};

    //for the Xi:
    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Xi0(1530)
    //3 = from Omega
    //4 = misidentified
    const unsigned NumChannels_Xim = 5;
    const double XimPurity = 0.9;
    //ratio Xi-(1530) to Xi0 (m=minus)
    const double Xim1530_to_Xim = 0.32*(1./3.);
    //ratio Xi0(1530) to Xi0 (n=neutral)
    const double Xin1530_to_Xim = 0.32*(2./3.);
    const double Omegam_to_Xim = 0.1;//ratio Omega/Xi
    const double OmegamXim_BR = 0.086;//branching ratio Omega->Xi
    const double Purities_Xim[NumChannels_Xim] = {XimPurity,XimPurity,XimPurity,XimPurity,1.-XimPurity};
    const double Fraction_Xim[NumChannels_Xim] = {1.-Xim1530_to_Xim-Xin1530_to_Xim-Omegam_to_Xim*OmegamXim_BR,Xim1530_to_Xim,Xin1530_to_Xim,Omegam_to_Xim*OmegamXim_BR,1};


    printf("Legend:\n");
    printf(" A : primary particle\n");
    printf(" A_B : feed-down from B to A\n");
    printf(" A' : missidentified particle\n");
    printf("\n");

    TString* ChannelName_p = new TString [NumChannels_p];
    ChannelName_p[0] = "p";
    ChannelName_p[1] = "p_Lambda";
    ChannelName_p[2] = "p_Sigma+";
    ChannelName_p[3] = "p'";

    TString* ChannelName_L = new TString [NumChannels_L];
    ChannelName_L[0] = "Lambda";
    ChannelName_L[1] = "Lambda_Sigma0";
    ChannelName_L[2] = "Lambda_Xi-";
    ChannelName_L[3] = "Lambda_Xi0";
    ChannelName_L[4] = "Lambda'";

    TString* ChannelName_Xim = new TString [NumChannels_Xim];
    ChannelName_Xim[0] = "Xi-";
    ChannelName_Xim[1] = "Xi-_Xi-(1530)";
    ChannelName_Xim[2] = "Xi-_Xi0(1530)";
    ChannelName_Xim[3] = "Xi-_Omega-";
    ChannelName_Xim[4] = "Xi-'";

    char* buffer = new char[128];

    double LambdaPar;
    double TotProb=0;
    printf("p_p:\n");
    for(unsigned uCh_p=0; uCh_p<NumChannels_p; uCh_p++){
        for(unsigned uCh_p2=0; uCh_p2<NumChannels_p; uCh_p2++){
            if(Compact && uCh_p>uCh_p2) continue;
            LambdaPar = Purities_p[uCh_p]*Fraction_p[uCh_p]*Purities_p[uCh_p2]*Fraction_p[uCh_p2];
            if(Compact && uCh_p!=uCh_p2) LambdaPar *= 2;
            TotProb+=LambdaPar;
            printf(" lambda(%s,%s) = %.2f%\n", ChannelName_p[uCh_p].Data(), ChannelName_p[uCh_p2].Data(), LambdaPar*100);
        }

    }
    //printf(" TotProb = %.3f\n", TotProb*100);
    printf("\np_Lambda:\n");
    TotProb=0;
    for(unsigned uCh_p=0; uCh_p<NumChannels_p; uCh_p++){
        for(unsigned uCh_L=0; uCh_L<NumChannels_L; uCh_L++){
            LambdaPar = Purities_p[uCh_p]*Fraction_p[uCh_p]*Purities_L[uCh_L]*Fraction_L[uCh_L];
            TotProb+=LambdaPar;
            printf(" lambda(%s,%s) = %.2f%\n", ChannelName_p[uCh_p].Data(), ChannelName_L[uCh_L].Data(), LambdaPar*100);
        }

    }
    //printf(" TotProb = %.3f\n", TotProb*100);
/*
    printf("SMALL TEST:\n");
    printf(" %f vs %f\n",   Purities_p[3]*Purities_L[0]+
                            Purities_p[0]*Purities_L[4]+
                            Purities_p[3]*Purities_L[4]
                            ,
                            Purities_p[3]+Purities_L[4]-
                            Purities_p[3]*Purities_L[4]
           );
*/
    printf("\nLambda_Lambda:\n");
    TotProb=0;
    for(unsigned uCh_L=0; uCh_L<NumChannels_L; uCh_L++){
        for(unsigned uCh_L2=0; uCh_L2<NumChannels_L; uCh_L2++){
            if(Compact && uCh_L>uCh_L2) continue;
            LambdaPar = Purities_L[uCh_L]*Fraction_L[uCh_L]*Purities_L[uCh_L2]*Fraction_L[uCh_L2];
            if(Compact && uCh_L!=uCh_L2) LambdaPar *= 2;
            TotProb+=LambdaPar;
            printf(" lambda(%s,%s) = %.2f%\n", ChannelName_L[uCh_L].Data(), ChannelName_L[uCh_L2].Data(), LambdaPar*100);
        }

    }
    //printf(" TotProb = %.3f\n", TotProb*100);

    printf("\np_Xi-:\n");
    TotProb=0;
    for(unsigned uCh_p=0; uCh_p<NumChannels_p; uCh_p++){
        for(unsigned uCh_Xim=0; uCh_Xim<NumChannels_Xim; uCh_Xim++){
            LambdaPar = Purities_p[uCh_p]*Fraction_p[uCh_p]*Purities_Xim[uCh_Xim]*Fraction_Xim[uCh_Xim];
            TotProb+=LambdaPar;
            printf(" lambda(%s,%s) = %.2f%\n", ChannelName_p[uCh_p].Data(), ChannelName_Xim[uCh_Xim].Data(), LambdaPar*100);
        }

    }

    delete [] ChannelName_p;
    delete [] ChannelName_L;
    delete [] ChannelName_Xim;
    delete [] buffer;

}
