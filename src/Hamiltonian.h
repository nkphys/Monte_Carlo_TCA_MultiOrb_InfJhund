#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void zheev_(char *, char *, int *, std::complex<double> *, int *, double *,
                       std::complex<double> *, int *, double *, int *);

class Hamiltonian
{
public:
    Hamiltonian(Parameters &Parameters__, Coordinates &Coordinates__, Coordinates &CoordinatesCluster__, MFParams &MFParams__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), CoordinatesCluster_(CoordinatesCluster__), MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
        HTBClusterCreate();
    }

    void Initialize();                                     //::DONE
    void Hoppings();                                       //::DONE
    double GetCLEnergy();                                  //::DONE
    void InteractionsCreate();                             //::DONE
    void InteractionsClusterCreate(int Center_site);       //::DONE
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE
    void HTBClusterCreate();                               //::DONE
    double chemicalpotential(double muin, double filling); //::DONE

    double chemicalpotentialCluster(double muin, double filling); //::DONE

    double TotalDensity();   //::DONE
    double ClusterDensity(); //::DONE
    double E_QM();           //::DONE

    double E_QMCluster();                 //::DONE
    void Diagonalize(char option);        //::DONE
    void DiagonalizeCluster(char option); //::DONE
    void copy_eigs(int i);                //::DONE
    void copy_eigs_Cluster(int i);        //::DONE

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    Coordinates &CoordinatesCluster_;
    MFParams &MFParams_;
    int lx_, ly_, ncells_, n_orbs_;
    int lx_cluster_, ly_cluster_, ncells_cluster;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> HTBCluster_;
    Matrix<complex<double>> Ham_;
    Matrix<complex<double>> HamCluster_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    vector<double> eigs_, eigsCluster_, eigsCluster_saved_, eigs_saved_, sx_, sy_, sz_;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp ;


    double HS_factor;
};

double Hamiltonian::chemicalpotential(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.005 * (eigs_[nstate - 1] - eigs_[0]) / nstate;
    N = filling * double(eigs_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (!Parameters_.fix_mu)
    {
        assert(!Parameters_.fix_mu);

        if (1 == 1)
        {
            for (int i = 0; i < 1000000; i++)
            {
               // cout<<mu_out<<"  "<<n1<<endl;
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigs_[j] - mu_out) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    mu_out += (N - n1) * dMubydN;
                    //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
                }
            }

            if (!converged)
            {
                cout << "mu_not_converged and N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }
        }

        double mu1, mu2;
        double mu_temp = muin;
        //cout<<"mu_input = "<<mu_temp<<endl;
        if (1 == 2)
        {
            mu1 = eigs_[0];
            mu2 = eigs_[nstate - 1];
            for (int i = 0; i < 4000000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigs_[j] - mu_temp) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    if (n1 > N)
                    {
                        mu2 = mu_temp;
                        mu_temp = 0.5 * (mu1 + mu_temp);
                    }
                    else
                    {
                        mu1 = mu_temp;
                        mu_temp = 0.5 * (mu2 + mu_temp);
                    }
                }
                //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

        return mu_out;
    }
    else
    {
        assert(Parameters_.fix_mu);
        return Parameters_.fixed_mu_value;
    }
} // ----------

double Hamiltonian::chemicalpotentialCluster(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigsCluster_.size();
    dMubydN = 0.001 * (eigsCluster_[nstate - 1] - eigsCluster_[0]) / nstate;
    N = filling * double(eigsCluster_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (!Parameters_.fix_mu)
    {
        assert(!Parameters_.fix_mu);
        if (1 == 1)
        {
            for (int i = 0; i < 100000; i++)
            {
                //cout<<mu_out<<"  "<<n1<<endl;
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigsCluster_[j] - mu_out) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    mu_out += (N - n1) * dMubydN;
                    //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
                }
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }
        }

        double mu1, mu2;
        double mu_temp = muin;
        //cout<<"mu_input = "<<mu_temp<<endl;
        if (1 == 2)
        {
            mu1 = eigsCluster_[0];
            mu2 = eigsCluster_[nstate - 1];
            for (int i = 0; i < 4000000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigsCluster_[j] - mu_temp) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    if (n1 > N)
                    {
                        mu2 = mu_temp;
                        mu_temp = 0.5 * (mu1 + mu_temp);
                    }
                    else
                    {
                        mu1 = mu_temp;
                        mu_temp = 0.5 * (mu2 + mu_temp);
                    }
                }
                //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

        return mu_out;
    }
    else
    {
        assert(Parameters_.fix_mu);
        return Parameters_.fixed_mu_value;
    }
} // ----------

void Hamiltonian::Initialize()
{

    ly_cluster_ = Parameters_.ly_cluster;
    lx_cluster_ = Parameters_.lx_cluster;
    ncells_cluster = lx_cluster_*ly_cluster_;

    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_*ly_;
    n_orbs_ = Parameters_.n_orbs;
    int space = ncells_ * n_orbs_;
    int spaceCluster = ncells_cluster* n_orbs_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
    HTBCluster_.resize(spaceCluster, spaceCluster);
    HamCluster_.resize(spaceCluster, spaceCluster);
    eigs_.resize(space);
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);
    eigs_saved_.resize(space);
    eigsCluster_.resize(spaceCluster);
    eigsCluster_saved_.resize(spaceCluster);

} // ----------

double Hamiltonian::TotalDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        n1 += 1.0 / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus * 1.0)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian::ClusterDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        n1 += 1.0 / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster * 1.0)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian::E_QM()
{
    double E = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        E += (eigs_[j]) / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian::E_QMCluster()
{
    double E = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        E += (eigsCluster_[j]) / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian::GetCLEnergy()
{

    double EClassical;
    int cell;
    double ei, ai;

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            cell = Coordinates_.Ncell(i, j);
            ei = MFParams_.etheta(i, j);
            ai = MFParams_.ephi(i, j);
            sx_[cell] = MFParams_.Moment_Size(i, j) * cos(ai) * sin(ei);
            sy_[cell] = MFParams_.Moment_Size(i, j) * sin(ai) * sin(ei);
            sz_[cell] = MFParams_.Moment_Size(i, j) * cos(ei);
        }
    }

    // Classical Energy
    EClassical = double(0.0);

    int _ix, _iy;
    for (int i = 0; i < ncells_; i++)
    {
        _ix = Coordinates_.indx_cellwise(i);
        _iy = Coordinates_.indy_cellwise(i);

        cell = Coordinates_.neigh(i, 0); //+x
        EClassical += 1.0 * Parameters_.K1x * ( (sx_[i] * sx_[cell]) + (sy_[i] * sy_[cell]) + (1.0 * sz_[i] * sz_[cell]));
        cell = Coordinates_.neigh(i, 2); //+y
        EClassical += Parameters_.K1y * ((sx_[i] * sx_[cell]) + (sy_[i] * sy_[cell]) + (1.0 * sz_[i] * sz_[cell]));
    }

    return EClassical;
} // ----------

void Hamiltonian::InteractionsCreate()
{


    //Kinetic Energy term
    int mx = Parameters_.TBC_mx;
    int my = Parameters_.TBC_my;

    complex<double> phasex, phasey;
    int l, m, a, b;
    double theta_l,phi_l, theta_m, phi_m;
    int lx_pos, ly_pos;
    int mx_pos, my_pos;

    HTB_.fill(0.0);

    for (l = 0; l < ncells_; l++)
    {
        lx_pos = Coordinates_.indx_cellwise(l);
        ly_pos = Coordinates_.indy_cellwise(l);
        theta_l = MFParams_.etheta(lx_pos, ly_pos);
        phi_l = MFParams_.ephi(lx_pos, ly_pos);


        // * +x direction Neighbor
        if (lx_pos == (Coordinates_.lx_ - 1))
        {
            phasex = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * mx) * PI / (1.0 * Parameters_.TBC_cellsX));
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 0); //+x neighbour cell
        mx_pos = Coordinates_.indx_cellwise(m);
        my_pos = Coordinates_.indy_cellwise(m);
        theta_m = MFParams_.etheta(mx_pos, my_pos);
        phi_m = MFParams_.ephi(mx_pos, my_pos);

        for(int orb1=0;orb1<n_orbs_;orb1++){
            for(int orb2=0;orb2<n_orbs_;orb2++){
                if(Parameters_.hopping_NN_X(orb1,orb2)!=0.0){
                    a = Coordinates_.Nbasis(lx_pos,ly_pos,orb1);
                    b = Coordinates_.Nbasis(mx_pos,my_pos,orb2);
                    assert(a != b);
                    if (a != b)
                    {
                        HTB_(b, a) = complex<double>(1.0 *Parameters_.hopping_NN_X(orb1,orb2), 0.0) * phasex *
                                ( (cos(theta_l*0.5)*cos(theta_m*0.5))  +  (sin(theta_l*0.5)*sin(theta_m*0.5)*exp(iota_complex*(phi_l-phi_m)))   );
                        HTB_(a, b) = conj(HTB_(b, a));
                    }
                }
            }
        }


        // * +y direction Neighbor
        if (ly_pos == (Coordinates_.ly_ - 1))
        {
            phasex = one_complex;
            phasey = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * my) * PI / (1.0 * Parameters_.TBC_cellsY));
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 2); //+y neighbour cell
        mx_pos = Coordinates_.indx_cellwise(m);
        my_pos = Coordinates_.indy_cellwise(m);
        theta_m = MFParams_.etheta(mx_pos, my_pos);
        phi_m = MFParams_.ephi(mx_pos, my_pos);

        for(int orb1=0;orb1<n_orbs_;orb1++){
            for(int orb2=0;orb2<n_orbs_;orb2++){
                if(Parameters_.hopping_NN_Y(orb1,orb2)!=0.0){

                    a = Coordinates_.Nbasis(lx_pos,ly_pos,orb1);
                    b = Coordinates_.Nbasis(mx_pos,my_pos,orb2);
                    assert(a != b);
                    if (a != b)
                    {
                        HTB_(b, a) = complex<double>(1.0*Parameters_.hopping_NN_Y(orb1,orb2), 0.0) * phasey *
                                ( (cos(theta_l*0.5)*cos(theta_m*0.5))  +  (sin(theta_l*0.5)*sin(theta_m*0.5)*exp(iota_complex*(phi_l-phi_m)))   );
                        HTB_(a, b) = conj(HTB_(b, a));
                    }
                }
            }
        }

    }




    //Interaction + onsite terms

    double ei, ai;
    int index;
    int i_posx, i_posy;

    Ham_ = HTB_;
    // Ham_.print();

    for (int i = 0; i < ncells_; i++)
    { // For each cell
        i_posx = Coordinates_.indx_cellwise(i);
        i_posy = Coordinates_.indy_cellwise(i);
        ei = MFParams_.etheta(i_posx, i_posy);
        ai = MFParams_.ephi(i_posx, i_posy);

        for(int orb=0;orb<n_orbs_;orb++){

            index=Coordinates_.Nbasis(i_posx, i_posy, orb);

            // On-Site potential
            a = Coordinates_.Nbasis(i_posx,i_posy,orb);
            Ham_(a, a) += complex<double>(1.0, 0.0) * (
                        Parameters_.OnSiteE[orb] +
                        MFParams_.Disorder(i_posx, i_posy)
                        );

        }
    }

} // ----------

void Hamiltonian::InteractionsClusterCreate(int Center_site)
{


    //Kinetic Energy term

    complex<double> phasex, phasey;
    int l, m, a, b;
    double theta_l,phi_l, theta_m, phi_m;
    int lx_pos, ly_pos;
    int mx_pos, my_pos;
    int x_pos, y_pos;
    int x_pos_p, y_pos_p;

    HTBCluster_.fill(0.0);

    for (l = 0; l < ncells_cluster; l++)
    {
        lx_pos = CoordinatesCluster_.indx_cellwise(l);
        ly_pos = CoordinatesCluster_.indy_cellwise(l);

        x_pos = Coordinates_.indx_cellwise(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx_cellwise(l);
        y_pos = Coordinates_.indy_cellwise(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy_cellwise(l);
        x_pos = (x_pos + Coordinates_.lx_) % Coordinates_.lx_;
        y_pos = (y_pos + Coordinates_.ly_) % Coordinates_.ly_;

        theta_l = MFParams_.etheta(x_pos, y_pos);
        phi_l = MFParams_.ephi(x_pos, y_pos);




        // * +x direction Neighbor
        m = CoordinatesCluster_.neigh(l, 0); //+x neighbour cell
        mx_pos = CoordinatesCluster_.indx_cellwise(m);
        my_pos = CoordinatesCluster_.indy_cellwise(m);
        x_pos_p = Coordinates_.indx_cellwise(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx_cellwise(m);
        y_pos_p = Coordinates_.indy_cellwise(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy_cellwise(m);
        x_pos_p = (x_pos_p + Coordinates_.lx_) % Coordinates_.lx_;
        y_pos_p = (y_pos_p + Coordinates_.ly_) % Coordinates_.ly_;

        theta_m = MFParams_.etheta(x_pos_p, y_pos_p);
        phi_m = MFParams_.ephi(x_pos_p, y_pos_p);


        if (mx_pos == (CoordinatesCluster_.lx_ - 1))
        {
            phasex = Parameters_.ClusterBoundaryConnection;
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }


        for(int orb1=0;orb1<n_orbs_;orb1++){
            for(int orb2=0;orb2<n_orbs_;orb2++){
                if(Parameters_.hopping_NN_X(orb1,orb2)!=0.0){
                    a = CoordinatesCluster_.Nbasis(lx_pos,ly_pos,orb1);
                    b = CoordinatesCluster_.Nbasis(mx_pos,my_pos,orb2);
                    assert(a != b);
                    if (a != b)
                    {
                        HTBCluster_(b, a) = complex<double>(1.0 *Parameters_.hopping_NN_X(orb1,orb2), 0.0) * phasex *
                                ( (cos(theta_l*0.5)*cos(theta_m*0.5))  +  (sin(theta_l*0.5)*sin(theta_m*0.5)*exp(iota_complex*(phi_l-phi_m)))   );
                        HTBCluster_(a, b) = conj(HTBCluster_(b, a));
                    }
                }
            }
        }


        // * +y direction Neighbor
        m = CoordinatesCluster_.neigh(l, 2); //+y neighbour cell
        mx_pos = CoordinatesCluster_.indx_cellwise(m);
        my_pos = CoordinatesCluster_.indy_cellwise(m);

        x_pos_p = Coordinates_.indx_cellwise(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx_cellwise(m);
        y_pos_p = Coordinates_.indy_cellwise(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy_cellwise(m);
        x_pos_p = (x_pos_p + Coordinates_.lx_) % Coordinates_.lx_;
        y_pos_p = (y_pos_p + Coordinates_.ly_) % Coordinates_.ly_;

        theta_m = MFParams_.etheta(x_pos_p, y_pos_p);
        phi_m = MFParams_.ephi(x_pos_p, y_pos_p);

        if (my_pos == (CoordinatesCluster_.ly_ - 1))
        {
            phasex = one_complex;
            phasey = Parameters_.ClusterBoundaryConnection;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }


        for(int orb1=0;orb1<n_orbs_;orb1++){
            for(int orb2=0;orb2<n_orbs_;orb2++){
                if(Parameters_.hopping_NN_Y(orb1,orb2)!=0.0){

                    a = CoordinatesCluster_.Nbasis(lx_pos,ly_pos,orb1);
                    b = CoordinatesCluster_.Nbasis(mx_pos,my_pos,orb2);
                    assert(a != b);
                    if (a != b)
                    {
                        HTBCluster_(b, a) = complex<double>(1.0*Parameters_.hopping_NN_Y(orb1,orb2), 0.0) * phasey *
                                ( (cos(theta_l*0.5)*cos(theta_m*0.5))  +  (sin(theta_l*0.5)*sin(theta_m*0.5)*exp(iota_complex*(phi_l-phi_m)))   );
                        HTBCluster_(a, b) = conj(HTBCluster_(b, a));
                    }
                }
            }
        }

    }

    //-----------------------------------


    int i_posx, i_posy;

    HamCluster_ = HTBCluster_;

    for (int i = 0; i < ncells_cluster; i++)
    { // For each cell in cluster

        i_posx = CoordinatesCluster_.indx_cellwise(i);
        i_posy = CoordinatesCluster_.indy_cellwise(i);

        x_pos = Coordinates_.indx_cellwise(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx_cellwise(i);
        y_pos = Coordinates_.indy_cellwise(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy_cellwise(i);
        x_pos = (x_pos + Coordinates_.lx_) % Coordinates_.lx_;
        y_pos = (y_pos + Coordinates_.ly_) % Coordinates_.ly_;

        for(int orb=0;orb<n_orbs_;orb++){

            a = CoordinatesCluster_.Nbasis(i_posx,i_posy,orb);
            HamCluster_(a, a) += complex<double>(1.0, 0.0) * (
                        Parameters_.OnSiteE[orb] +
                        MFParams_.Disorder(x_pos, y_pos)
                        );

        }
    }


} // ----------

void Hamiltonian::Check_up_down_symmetry()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < ncells_*n_orbs_; i++)
    {
        for (int j = 0; j < ncells_*n_orbs_; j++)
        {
            temp2 = Ham_(i, j) - Ham_(i + ncells_*n_orbs_, j + ncells_*n_orbs_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp += temp2 * conj(temp2);
        }
    }

    cout << "Assymetry in up-down sector: " << temp << endl;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < HamCluster_.n_row(); i++)
    {
        for (int j = 0; j < HamCluster_.n_row(); j++)
        {
            if (HamCluster_(i, j) != conj(HamCluster_(j, i)))
            {
                cout << i << "," << j << endl;
                cout << "i,j = " << HamCluster_(i, j) << ", j,i=" << conj(HamCluster_(j, i)) << endl;
            }
            assert(HamCluster_(i, j) == conj(HamCluster_(j, i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}

void Hamiltonian::Diagonalize(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = Ham_.n_row();
    int lda = Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(), eigs_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian::DiagonalizeCluster(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    // cout << option;
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = HamCluster_.n_row();
    int lda = HamCluster_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigsCluster_.resize(HamCluster_.n_row());
    fill(eigsCluster_.begin(), eigsCluster_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian::HTBCreate()
{


    //Convention used
    //orb=0=d
    //orb=1=px
    //orb=2=py





} // ----------

void Hamiltonian::HTBClusterCreate()
{


    // HTBCluster_.print();

} // ----------

void Hamiltonian::Hoppings()
{
    //Using matrices from Parameters_

} // ----------

void Hamiltonian::copy_eigs(int i)
{

    int space =  ncells_ *n_orbs_;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigs_saved_[j] = eigs_[j];
        }
    }
}

void Hamiltonian::copy_eigs_Cluster(int i)
{

    int space = ncells_cluster *n_orbs_;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_[j] = eigsCluster_saved_[j];
        }
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_saved_[j] = eigsCluster_[j];
        }
    }
}

#endif
