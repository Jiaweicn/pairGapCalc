//using mass to calc gap
//for proton gap, neutron gap and neutron-proton interaction energy
/*see paper doi.org/10.1016/0375-9474(88)90370-3*/

#include "../constant.h"
#include "ameMass2020.h"
#include "shared.h"

void initMGR(TGraph *grs[3][2],TMultiGraph *mgr){
	const short markers[]={71,73};//open
	const short colors[]={4,2,1};
	const short lines[]={2,5};
	for(short i=0;i<3;i++){//different kinds of nucleon-gaps
		for(short j=0;j<2;j++){//different cantegories of nuclei
			grs[i][j]=new TGraph();//different cantegories of nuclei
			grs[i][j]->SetMarkerSize(2);
			grs[i][j]->SetLineWidth(2);
			grs[i][j]->SetMarkerStyle(markers[j]);
			grs[i][j]->SetMarkerColor(colors[i]);
			grs[i][j]->SetLineColor(colors[i]);
			grs[i][j]->SetLineStyle(lines[j]);
			if(i==0){
				grs[i][j]->SetTitle("neutron pair gap");
				grs[i][j]->GetXaxis()->SetTitle("Neutron number N");
				grs[i][j]->GetYaxis()->SetTitle("#Delta_{nn}/MeV");
			}
			else if(i==1){
				grs[i][j]->SetTitle("proton pair gap");
				grs[i][j]->GetXaxis()->SetTitle("Proton number Z");
				grs[i][j]->GetYaxis()->SetTitle("#Delta_{pp}/MeV");
			}
			else if(i==2){
				grs[i][j]->SetTitle("neutron-proton interaction energy");
				grs[i][j]->GetXaxis()->SetTitle("Mass number A");
				grs[i][j]->GetYaxis()->SetTitle("#delta_{np}/MeV");
			}
			mgr->Add(grs[i][j]);
		}
	}
}
void pairGap(short z0=50,short n0=0){//3 kind of nucleon-pair gaps of a series of isotopes/isotones
	short a0=0;
	if(z0>0 && n0==0){//isotopes
		a0=z0;
		cout<<"calculating gaps for isotopes: Z="<<z0<<"\n";
	}
	else if(n0>0 && z0==0){//isotones
		cout<<"calculating gaps for isotones: N="<<n0<<"\n";
		a0=n0;
	}
	else{
		cout<<"wrong input!\n";
		exit(1);
	}
	double ppGap,nnGap,npDelta;
	TGraph *grs[3][2];
	auto *mgr=new TMultiGraph();
	initMGR(grs,mgr);
	short type;
	string elementName;
	bool nameGot=false;

	for(short a=a0;a<300;a++){//mass number
		if(z0>0 && n0==0){//isotopes
			const short n=a-z0;
			if(!nucleiType(z0,n,type)) continue;
			//if(calcNPgap1993(a,z0,npDelta))//np interaction energy, MeV
			if(calcNPgap2020(a,z0,npDelta))//np interaction energy, MeV
				grs[2][type%2]->SetPoint(grs[2][type%2]->GetN(),n,npDelta);
			else continue;
			//if(calcNgap1993(a,z0,nnGap,npDelta))//nn pair gap, MeV
			if(calcNgap2020(a,z0,nnGap,npDelta))//nn pair gap, MeV
				grs[0][type%2]->SetPoint(grs[0][type%2]->GetN(),n,nnGap);
			//if(calcPgap1993(a,z0,ppGap,npDelta))//pp pair gap, MeV
			if(calcPgap2020(a,z0,ppGap,npDelta))//pp pair gap, MeV
				grs[1][type%2]->SetPoint(grs[1][type%2]->GetN(),n,ppGap);
			if(!nameGot){
				nameGot=getNucleusName1993(a,z0,elementName);
				nameGot=true;
			}
		}
		if(n0>0 && z0==0){//isotones
			const short z=a-n0;
			if(!nucleiType(z,n0,type)) continue;
			//if(calcNPgap1993(a,z,npDelta))//np interaction energy, MeV
			if(calcNPgap2020(a,z,npDelta))//np interaction energy, MeV
				grs[2][type/2]->SetPoint(grs[2][type/2]->GetN(),z,npDelta);
			//if(calcNgap1993(a,z,nnGap,npDelta))//nn pair gap, MeV
			if(calcNgap2020(a,z,nnGap,npDelta))//nn pair gap, MeV
				grs[0][type/2]->SetPoint(grs[0][type/2]->GetN(),z,nnGap);
			//if(calcPgap1993(a,z,ppGap,npDelta))//pp pair gap, MeV
			if(calcPgap2020(a,z,ppGap,npDelta))//pp pair gap, MeV
				grs[1][type/2]->SetPoint(grs[1][type/2]->GetN(),z,ppGap);
		}
	}
	auto *c=new TCanvas("cgap","cgap",1200,800);
	auto *lgde=new TLegend(0.65,0.75,0.8,0.9);
	lgde->SetFillStyle(0);
	auto *lgdo=new TLegend(0.8,0.75,0.9,0.9);
	lgdo->SetFillStyle(0);
	if(z0>0 && n0==0){
		lgde->AddEntry(grs[0][0],"#Delta_{nn}, eN");
		lgde->AddEntry(grs[1][0],"#Delta_{pp}, eN");
		lgde->AddEntry(grs[2][0],"#delta_{np}, eN");
		lgdo->AddEntry(grs[0][1],"oN");
		lgdo->AddEntry(grs[1][1],"oN");
		lgdo->AddEntry(grs[2][1],"oN");
	}
	else if(n0>0 && z0==0){
		lgde->AddEntry(grs[0][0],"#Delta_{nn}, eZ");
		lgde->AddEntry(grs[1][0],"#Delta_{pp}, eZ");
		lgde->AddEntry(grs[2][0],"#delta_{np}, eZ");
		lgdo->AddEntry(grs[0][1],"oZ");
		lgdo->AddEntry(grs[1][1],"oZ");
		lgdo->AddEntry(grs[2][1],"oZ");
	}

	mgr->Draw("ap");
	if(z0>0 && n0==0){
		mgr->SetTitle(Form("nucleon-pair gaps of %s(Z=%d), #Delta(#delta)'s",elementName.c_str(),z0));
		mgr->GetXaxis()->SetTitle("Neutron number N");
	}
	else if(n0>0 && z0==0){
		mgr->SetTitle(Form("nucleon-pair gaps of isotones(N=%d), #Delta(#delta)'s",n0));
		mgr->GetXaxis()->SetTitle("Proton number Z");
	}
	mgr->GetYaxis()->SetTitle("#Delta/MeV");
	lgde->Draw();
	lgdo->Draw();
}
