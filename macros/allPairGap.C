//calculation pair gap of ALL known isotopes
//using binding energy to calc gap
//using data from mass table of artemis
//only for proton gap and neutron gap
#include "/home/cai/prjs/include/cpp/outColor.h"
#include "../constant.h"
#include "ameMass2020.h"
#include "shared.h"

void initMGR(TGraph *grs[3][4],TMultiGraph *mgrs[3]){
	const short colors[]={1,2,4,6};
	const short markers[]={20,21,22,33};
	for(short i=0;i<3;i++){//different kinds of nucleon-gaps
		mgrs[i]=new TMultiGraph();
		for(short j=0;j<4;j++){//different cantegories of nuclei
			grs[i][j]=new TGraph();//different cantegories of nuclei
			grs[i][j]->SetMarkerSize(1);
			grs[i][j]->SetMarkerStyle(markers[j]);
			grs[i][j]->SetMarkerColor(colors[j]);
			grs[i][j]->SetLineColor(colors[j]);
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
			mgrs[i]->Add(grs[i][j]);
		}
	}
}
void pairGap(const short &z0,TGraph *grs[3][4]){//pair gap of the isotopes of an element
	double ppGap,nnGap,npDelta;
	short type;
	const float xshift=0.;//shift of x
	for(short a=z0;a<300;a++){//mass number starts from charge number
		const short n=a-z0;
		if(!nucleiType(z0,n,type)) continue;
		//if(calcNPgap1993(a,z0,npDelta))//np interaction energy, MeV
		if(calcNPgap2020(a,z0,npDelta))//np interaction energy, MeV
			grs[2][type]->SetPoint(grs[2][type]->GetN(),a+xshift*type,npDelta);
		else continue;
		//if(calcNgap1993(a,z0,nnGap,npDelta))//nn pair gap, MeV
		if(calcNgap2020(a,z0,nnGap,npDelta))//nn pair gap, MeV
			grs[0][type]->SetPoint(grs[0][type]->GetN(),n+xshift*type,nnGap);
		//if(calcPgap1993(a,z0,ppGap,npDelta))//pp pair gap, MeV
		if(calcPgap2020(a,z0,ppGap,npDelta))//pp pair gap, MeV
			grs[1][type]->SetPoint(grs[1][type]->GetN(),z0+xshift*type,ppGap);
	}
}

void allPairGap(){//nucleon-pair gaps
	TGraph *grs[3][4];
	TMultiGraph *mgrs[3];
	initMGR(grs,mgrs);
	auto *lgd=new TLegend(0.78,0.7,0.9,0.9);
	lgd->SetFillStyle(0);
	lgd->AddEntry(grs[0][0],"eZ-eN");
	lgd->AddEntry(grs[0][1],"eZ-oN");
	lgd->AddEntry(grs[0][2],"oZ-eN");
	lgd->AddEntry(grs[0][3],"oZ-oN");

	for(short z=1;z<120;z++) pairGap(z,grs);

	const float ylows[]={-0.1,0.4,-1.1};
	const float yhighs[]={5,5,3};
	TCanvas *cgaps[3];
	for(short i=0;i<3;i++){
		cgaps[i]=new TCanvas("","",900,600);
		mgrs[i]->Draw("ap");
		/*if(i==0) mgrs[i]->GetXaxis()->SetLimits(0,160);
		else if(i==1) mgrs[i]->GetXaxis()->SetLimits(0,120);
		else mgrs[i]->GetXaxis()->SetLimits(0,280);*/
		mgrs[i]->GetHistogram()->GetYaxis()->SetRangeUser(ylows[i],yhighs[i]);
		mgrs[i]->GetHistogram()->SetTitle(grs[i][0]->GetTitle());
		mgrs[i]->GetXaxis()->SetTitle(grs[i][0]->GetXaxis()->GetTitle());
		mgrs[i]->GetYaxis()->SetTitle(grs[i][0]->GetYaxis()->GetTitle());
		lgd->Draw();
	}
}
