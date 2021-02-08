#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams
{
public:
    // Define Fields
    Matrix<double> etheta, ephi;
    Matrix<double> Sz, Sx, Sy;
    Matrix<double> etheta_avg, ephi_avg;
    Matrix<double> Moment_Size;
    Matrix<double> Disorder;

    // Constructor
    MFParams(Parameters &Parameters__, Coordinates &Coordinates__, mt19937_64 &Generator1__, mt19937_64 &Generator2__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }

    double random1();
    double random2();
    void FieldThrow(int site, string mc_dof_type);
    void initialize();
    void Adjust_MCWindow();
    void Calculate_Fields_Avg();
    void Read_classical_DOFs(string filename);
    void Push_to_Prob_Distributions(double theta_, double phi_);
    double Lorentzian(double x, double brd);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int lx_, ly_, ns_;

    uniform_real_distribution<double> dis1_; //for random fields
    uniform_real_distribution<double> dis2_; //for random disorder

    Mat_1_doub Distribution_Phi, Distribution_Theta;
    double d_Phi, d_Theta;
    int N_Phi, N_Theta;

    //mt19937_64 mt_rand(Parameters_.RandomSeed);
};

void MFParams::Adjust_MCWindow()
{
    double ratio;
    ratio = Parameters_.AccCount[0] / (Parameters_.AccCount[0] + Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0] = 0;
    Parameters_.AccCount[1] = 0;
    Parameters_.WindowSize *= abs(1.0 + 0.1 * (ratio - 0.5));
    if(Parameters_.WindowSize>1){
        Parameters_.WindowSize=1.0;
    }
    // Parameters_.WindowSize =0.2;
    cout << "Ratio: " << ratio << "  window size:  " << Parameters_.WindowSize << endl;
    return;
} // ----------

void MFParams::FieldThrow(int site, string mc_dof_type)
{
    int a, b;

    int Pi_multiple;

    double Pi = Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    a = Coordinates_.indx_cellwise(site);
    b = Coordinates_.indy_cellwise(site);

    //ANGLES
    if (mc_dof_type == "phi")
    {
        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;

        // Pi_multiple = ephi(a, b)/Pi;


        //        if (ephi(a, b) < 0.0)
        //        {
        //            ephi(a, b) = -ephi(a, b);
        //        }

        //
        ephi(a, b) = fmod(ephi(a, b), 2.0 * Pi);
        //-------------
        if( ephi(a,b) < 0.0) {ephi(a,b) += 2.0*Pi; }
        if( ephi(a,b) >=2.0*Pi) {ephi(a,b) -= 2.0*Pi;}
        //-----------
    }

    if (mc_dof_type == "theta")
    {

        cout <<"Only theta is not allowed, use only phi instead"<<endl;
        assert(false);

    }

    if (mc_dof_type == "theta_ising")
    {

        if(abs(etheta(a, b))<0.0001){ //i.e.=0
            etheta(a, b) = Pi;
        }
        else{  //i.e. Pi
            if((abs(etheta(a, b)-Pi)>0.0001)){
                cout<<etheta(a, b)<<endl;
                assert(abs(etheta(a, b)-Pi)<0.0001);
            }

            etheta(a, b) = 0;
        }

    }


    if (mc_dof_type == "theta_and_phi")
    {


        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;
        if( ephi(a,b) < 0.0) {ephi(a,b) += 2.0*Pi; }
        if( ephi(a,b) >=2.0*Pi) {ephi(a,b) -= 2.0*Pi;}


        etheta(a, b) += Pi * (random1() - 0.5) * MC_Window;
        if ( etheta(a,b) < 0.0 ) {
            etheta(a,b) = - etheta(a,b);
            ephi(a,b) = fmod( ephi(a,b)+Pi, 2.0*Pi );
        }
        if ( etheta(a,b) > Pi ) {
            etheta(a,b) = 2.0*Pi - etheta(a,b);
            ephi(a,b) = fmod( ephi(a,b) + Pi, 2.0*Pi );
        }

    }


} // ----------

double MFParams::random1()
{

    return dis1_(Generator1_);
}

double MFParams::random2()
{

    return dis2_(Generator2_);
}

void MFParams::initialize()
{


    bool Diagonal_ZigZag_Ising_alongZ=false;
    bool two_by_two_Plaquettes_Ising_alongZ=false;
    bool FM_state_Ising=false;
    bool AFM_state_Ising=false;
    lx_ = Coordinates_.lx_;
    ly_ = Coordinates_.ly_;

    // srand(Parameters_.RandomSeed);


    etheta_avg.resize(lx_, ly_);
    ephi_avg.resize(lx_, ly_);
    Disorder.resize(lx_, ly_);

    etheta.resize(lx_, ly_);
    ephi.resize(lx_, ly_);
    Moment_Size.resize(lx_, ly_);

    Sz.resize(lx_, ly_);
    Sx.resize(lx_, ly_);
    Sy.resize(lx_, ly_);

    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file << "#seed=" << Parameters_.RandomDisorderSeed << " for mt19937_64 Generator is used" << endl;
    Disorder_conf_file << "#ix   iy    Dis[ix,iy]" << endl;

    ofstream Initial_MC_DOF_file("Initial_MC_DOF_values");

    Initial_MC_DOF_file << "#seed=" << Parameters_.RandomSeed << " for mt19937_64 Generator is used" << endl;
    Initial_MC_DOF_file << "#ix   iy   Theta(x,y)    Phi(x,y)      Moment_Size(x,y)" << endl;




    string temp_string;
    int ix_, iy_, spin_offset;
    if (Parameters_.Read_Seed_from_file_ == true)
    {
        ifstream Initial_Seed(Parameters_.Seed_file_name_);
        getline(Initial_Seed, temp_string);
        // cout << temp_string << endl;
        for (int ix = 0; ix < lx_; ix++)
        {
            for (int iy = 0; iy < ly_; iy++)
            {
                Initial_Seed >> ix_ >> iy_ >> etheta(ix, iy) >> ephi(ix, iy) >> Moment_Size(ix, iy);
                assert(ix_ == ix);
                assert(iy_ == iy);
            }
        }
    }

    else
    {

        //Initialization
        for (int j = 0; j < ly_; j++)
        {
            for (int i = 0; i < lx_; i++)
            {
                ephi(i, j) = 0.0;
                etheta(i, j) = PI*0.5;
            }
        }

        for (int j = 0; j < ly_; j++)
        {
            for (int i = 0; i < lx_; i++)
            {
                //RANDOM fields
                if (Parameters_.MC_on_theta_and_phi == true)
                {
                    ephi(i, j) = 2.0 * random1() * PI;
                    etheta(i, j) = random1() * PI;
                }
                else
                {
                    if (Parameters_.MC_on_phi == true)
                    {
                        ephi(i, j) = 2.0 * random1() * PI;
                    }


                    if (Parameters_.MC_on_theta == true)
                    {
                        if(Parameters_.ThetaIsing){
                            etheta(i, j) = int(random1()+0.5) * PI;
                        }
                        else{
                            etheta(i, j) = random1() * PI;
                        }
                    }


                    if( !Parameters_.MC_on_phi && !Parameters_.MC_on_theta){

                        if(Diagonal_ZigZag_Ising_alongZ){

                            if( ((i%4)==0) || ((i%4)==1)){
                                spin_offset=1;
                            }
                            else{
                                spin_offset=-1;
                            }


                            if(j%2==0){
                                iy_=j/2;
                            }
                            else{
                                iy_= (j -1)/2;
                            }


                            if( (i%4 == 0) ||  (i%4 == 2) ){


                                if(iy_%2==0){
                                    spin_offset = 1*spin_offset;
                                }
                                else{
                                    spin_offset = -1*spin_offset;
                                }

                            }
                            else{

                                if( (iy_%2 == 0) && (j%2==0) ){
                                    spin_offset = 1*spin_offset;
                                }
                                else if((iy_%2 == 1) && (j%2==0)){
                                    spin_offset = -1*spin_offset;
                                }
                                else if((iy_%2 == 0) && (j%2==1)){
                                    spin_offset = -1*spin_offset;
                                }
                                else{
                                    assert ((iy_%2 == 1) && (j%2==1));
                                    spin_offset = 1*spin_offset;
                                }
                            }

                            etheta(i, j) = ((-1*spin_offset*1.0) + 1.0) *0.5* PI;
                        }


                        if(two_by_two_Plaquettes_Ising_alongZ){
                            if( ((i%4)==0) || ((i%4)==1)){
                                spin_offset=1;
                            }
                            else{
                                spin_offset=-1;
                            }


                            if(j%2==0){
                                iy_=j/2;
                            }
                            else{
                                iy_= (j -1)/2;
                            }


                            if(iy_%2==0){
                                spin_offset = 1*spin_offset;
                            }
                            else{
                                spin_offset = -1*spin_offset;
                            }



                            etheta(i, j) = ((-1*spin_offset*1.0) + 1.0) *0.5* PI;
                        }

                        if(FM_state_Ising){
                            etheta(i, j) = ((-1*1.0) + 1.0) *0.5* PI;
                        }
                        if(AFM_state_Ising){

                            spin_offset = int(pow(-1.0, i+j));
                            etheta(i, j) = ((-1*spin_offset*1.0) + 1.0) *0.5* PI;

                        }

                    }
                }


                Moment_Size(i, j) = 1.0;

                Initial_MC_DOF_file << i << setw(15) << j << setw(15) << etheta(i, j) << setw(15) << ephi(i, j)
                                    << setw(15) << Moment_Size(i, j) << endl;

            }
        }
    }

    //RANDOM Disorder
    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {
            Disorder(i, j) = Parameters_.Disorder_Strength * ((2.0 * random2()) - 1.0);
            Disorder_conf_file << i << "  " << j << "  " << Disorder(i, j) << endl;
        }
        Disorder_conf_file << endl;
    }



    //Initializing PDF's for Phi and Theta.
    N_Phi=2000;
    N_Theta=1000;
    d_Phi=(2.0*PI)/(1.0*N_Phi);
    d_Theta=PI/(1.0*N_Theta);
    Distribution_Phi.resize(N_Phi);
    Distribution_Theta.resize(N_Theta);
    //Phi_val=d_Phi*N_Phi
    //Theta_val=d_Theta*N_Theta;


} // ----------

double MFParams::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}

void MFParams::Push_to_Prob_Distributions(double theta_, double phi_){
    double eta_=0.01;

    for(int theta_no=0;theta_no<Distribution_Theta.size();theta_no++){
        Distribution_Theta[theta_no] +=Lorentzian(theta_-theta_no*d_Theta, eta_)*(1.0/(Parameters_.IterMax*Coordinates_.ncells_));
    }

    for(int phi_no=0;phi_no<Distribution_Phi.size();phi_no++){
        Distribution_Phi[phi_no] +=Lorentzian(phi_-phi_no*d_Phi, eta_)*(1.0/(Parameters_.IterMax*Coordinates_.ncells_));
    }

    //cout<<"theta = "<<theta_<<", phi = "<<phi_<<endl;

}

void MFParams::Calculate_Fields_Avg()
{

    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {

            ephi_avg(i, j) = ephi_avg(i, j) + ephi(i, j);
            etheta_avg(i, j) = etheta_avg(i, j) + etheta(i, j);
        }
    }

} // ----------

void MFParams::Read_classical_DOFs(string filename)
{

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    getline (fl_in,tmp_str);

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            fl_in >> tmp_double >> tmp_double >> etheta(i, j) >> ephi(i, j)>>tmp_double;
        }
    }


} // ----------

#endif
