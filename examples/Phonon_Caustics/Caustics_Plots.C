//This is a root Program to Plot the Caustic Patterns from the phonon_hits.txt File
//Function to Plot Only Phonon Caustic.
void TransFast(){
  ifstream in;
  TString B= "phonon_hits.txt";
  in.open("phonon_hits.txt");
   TH2D *Caustics= new TH2D("Caustics","Phonon Caustics",500,-0.002,0.002,500,-0.002,0.002);
  Int_t nlines = 0;
  TString Name_Phonon;
   Double_t X_f,Y_f,Z_f;


  while (1) {

     in>>Name_Phonon>>X_f>>Y_f>>Z_f;
     if(Name_Phonon=="phononTF"){
     Caustics->Fill(X_f,Y_f);
    }
     if (!in.good()) break;
  nlines++;
  }

  TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
  c1->SetFillColor(1);
  Caustics->SetMarkerColorAlpha(kWhite, 0.2);
  //c1->SetLogz();

  //Caustics->Draw("colz");
  Caustics->Draw();
  in.close();

}

void TransSlow(){
  ifstream in;
  TString B= "phonon_hits.txt";
  in.open("phonon_hits.txt");
   TH2D *Caustics= new TH2D("Caustics","Phonon Caustics",500,-0.002,0.002,500,-0.002,0.002);
  Int_t nlines = 0;
  TString Name_Phonon;
   Double_t X_f,Y_f,Z_f;


  while (1) {

     in>>Name_Phonon>>X_f>>Y_f>>Z_f;
     if(Name_Phonon=="phononTS"){
     Caustics->Fill(X_f,Y_f);
    }
     if (!in.good()) break;
  nlines++;
  }
  //printf(" found %d points\n",nlines);
  //gStyle->SetPalette(kBlack+3);
  TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
  c1->SetFillColor(1);
  Caustics->SetMarkerColorAlpha(kWhite, 0.2);
  //c1->SetLogz();

  //Caustics->Draw("colz"); If you want to see on a Blue and yellow colors
  Caustics->Draw();
  in.close();

}

void TransFast_and_Slow(){
  ifstream in;
  TString B= "phonon_hits.txt";
  in.open("phonon_hits.txt");
   TH2D *Caustics= new TH2D("Caustics","Phonon Caustics",500,-0.002,0.002,500,-0.002,0.002);
  Int_t nlines = 0;
  TString Name_Phonon;
   Double_t X_f,Y_f,Z_f;


  while (1) {

     in>>Name_Phonon>>X_f>>Y_f>>Z_f;

     Caustics->Fill(X_f,Y_f);

     if (!in.good()) break;
  nlines++;
  }
  //printf(" found %d points\n",nlines);
  //gStyle->SetPalette(kBlack+3);
  TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
  //c1->SetFillColor(1);
  //Caustics->SetMarkerColorAlpha(kWhite, 0.2);
  //c1->SetLogz();

  Caustics->Draw("colz");
  //Caustics->Draw();
  in.close();

}



#include <iostream>
using namespace std;
void Caustics_Plots(TString Phonon_Name){

TString  Phonon_Case = Phonon_Name;
Int_t What_pnonon1,What_pnonon2,What_pnonon3;
What_pnonon1=Phonon_Case.CompareTo("Fast");
What_pnonon2=Phonon_Case.CompareTo("Slow");
What_pnonon3=Phonon_Case.CompareTo("Both");
if (What_pnonon1==0) {
cout<<"Fast";
TransFast();

}
else if (What_pnonon2==0){cout<<"Slow";TransSlow(); }
else {cout<<"Both";TransFast_and_Slow();}

}
