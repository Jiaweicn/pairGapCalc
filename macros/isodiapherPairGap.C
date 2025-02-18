#include "../constant.h"
#include "ameMass2020.h"
#include "shared.h"

void initMGR(const short nMinusZ,TGraph *grs[3][4],TMultiGraph *mgrs[3]){
	const short colors[]={4,2,1,6};
	const short markers[]={89,90,91,95};
	const short lines[]={9,5,7,2};
	for(short i=0;i<3;i++){//different cantegories of nuclei
		mgrs[i]=new TMultiGraph();
		for(short j=0;j<4;j++){//different cantegories of nuclei
			grs[i][j]=new TGraph();
			grs[i][j]->SetMarkerSize(3);
			grs[i][j]->SetLineWidth(3);
			grs[i][j]->SetMarkerStyle(markers[j]);
			grs[i][j]->SetLineStyle(lines[j]);
			grs[i][j]->SetMarkerColor(colors[i]);
			grs[i][j]->SetLineColor(colors[i]);
			if(i==0){//nn
				grs[i][j]->SetTitle(Form("neutron-neutron pair gap(N-Z=%d)",nMinusZ));
				grs[i][j]->GetXaxis()->SetTitle("Neutron number N");
				grs[i][j]->GetYaxis()->SetTitle("#Delta_{nn}/MeV");
			}
			else if(i==1){//pp
				grs[i][j]->SetTitle(Form("proton-proton pair gap(N-Z=%d)",nMinusZ));
				grs[i][j]->GetXaxis()->SetTitle("Proton number Z");
				grs[i][j]->GetYaxis()->SetTitle("#Delta_{pp}/MeV");
			}
			else if(i==2){//np
				grs[i][j]->SetTitle(Form("neutron-proton interaction energy(N-Z=%d)",nMinusZ));
				grs[i][j]->GetXaxis()->SetTitle("Mass number A");
				//grs[i][j]->GetXaxis()->SetTitle("Proton number Z");
				grs[i][j]->GetYaxis()->SetTitle("#delta_{np}/MeV");
			}
			//mgrs[i]->Add(grs[i][j]);
		}
	}
}
void plot2(const short &nMinusZ,TGraph *grs[3][4]){
	string pairLgd[3]={"#Delta_{nn}","#Delta_{pp}","#delta_{np}"};
	string typeLgd[4]={"eZ-eN","eZ-oN","oZ-eN","oZ-oN"};
	auto *cgap=new TCanvas("cgap","cgap",900,600);
	auto *lgd=new TLegend(0.65,0.65,0.9,0.9);
	lgd->SetFillStyle(0);
	auto *mgr=new TMultiGraph();
	for(short i=0;i<2;i++){//different pairs
		for(short j=0;j<4;j++){//different cantegories of nuclei
			if(grs[i][j]->GetN()>0){
				mgr->Add(grs[i][j]);
				lgd->AddEntry(grs[i][j],(pairLgd[i]+": "+typeLgd[j]).c_str());
			}
		}
	}
	mgr->Draw("apl");
	mgr->SetTitle(Form("nucleon-pair gaps(N-Z=%d)",nMinusZ));
	mgr->GetXaxis()->SetTitle("p/n number");
	mgr->GetYaxis()->SetTitle("#Delta/MeV");
	lgd->Draw();
}
void plot3(const short &nMinusZ,TGraph *grs[3][4],TMultiGraph *mgrs[3]){
	string typeLgd[4]={"eZ-eN","eZ-oN","oZ-eN","oZ-oN"};
	TCanvas *cgaps[3];
	TLegend *lgds[3];
	for(short i=0;i<3;i++){//different pairs
		cgaps[i]=new TCanvas("","",900,600);
		lgds[i]=new TLegend(0.75,0.75,0.9,0.9);
		lgds[i]->SetFillStyle(0);
		for(short j=0;j<4;j++){//different cantegories of nuclei
			if(grs[i][j]->GetN()>0){
				mgrs[i]->Add(grs[i][j]);
				lgds[i]->AddEntry(grs[i][j],typeLgd[j].c_str());
			}
		}
		mgrs[i]->Draw("apl");
		mgrs[i]->SetTitle(grs[i][0]->GetTitle());
		mgrs[i]->GetXaxis()->SetTitle(grs[i][0]->GetXaxis()->GetTitle());
		mgrs[i]->GetYaxis()->SetTitle(grs[i][0]->GetYaxis()->GetTitle());
		lgds[i]->Draw();
	}
}
void isodiapherPairGap(const short nMinusZ=20){//pair gaps of an isodiapher chain
	TGraph *grs[3][4];
	TMultiGraph *mgrs[3];
	initMGR(nMinusZ,grs,mgrs);

	for(short z=1;z<120;z++){//charge of nuclei
		const short n=z+nMinusZ;//neutron number
		const short a=z+n;//mass number
		double nnGap,ppGap,npDelta;
		short type;
		if(!nucleiType(z,n,type)) continue;
		//if(calcNPgap1993(a,z,npDelta))//np pair gap, MeV
		if(calcNPgap2020(a,z,npDelta))//np pair gap, MeV
			grs[2][type]->SetPoint(grs[2][type]->GetN(),a,npDelta);
		else continue;
		//if(calcNgap1993(a,z,nnGap,npDelta))//nn pair gap, MeV
		if(calcNgap2020(a,z,nnGap,npDelta))//nn pair gap, MeV
			grs[0][type]->SetPoint(grs[0][type]->GetN(),n,nnGap);
		//if(calcPgap1993(a,z,ppGap,npDelta))//pp pair gap, MeV
		if(calcPgap2020(a,z,ppGap,npDelta))//pp pair gap, MeV
			grs[1][type]->SetPoint(grs[1][type]->GetN(),z,ppGap);
	}
	plot2(nMinusZ,grs);//plot nn & pp gaps in one canvas
	plot3(nMinusZ,grs,mgrs);//plot nn, pp & np gaps in 3 canvases
}
