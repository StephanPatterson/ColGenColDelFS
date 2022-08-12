//
//  main.cpp
//  ColGenFS
//
//  Created by StephanP on 9/18/20.
//  Copyright © 2020 Stephan Patterson. All rights reserved.
//
//  Using the General formulation directly and introducing (all, one, or some constant number of) columns h where ch - u'A is negative
//  This version is for deleting columns as they leave the basis

#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include "SuppPt.hpp"
#include <numeric>

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);
//std::vector< double> mvmult(const std::vector< std::vector<int> > &, const int, const long int, const long int, const std::vector< double> &);
//void mvmult(const std::vector< std::vector<int> > &, const int, const long int, const long int, const std::vector< double> &, double *);
void mvmult(const int , const int long , const int long , const std::vector< double> &, double * const, const int * const , const int );
void mvmult(const int , const int long , const int long , const std::vector< double> &, std::vector<int> &, double * const, const int * const , const int );
std::vector< double> mvmult(std::vector< std::vector<int> >::iterator,const int, const std::vector< double> &);
std::vector< double> vmmult(const std::vector< std::vector<int> > & , const int , const long int , const double* const );
int vmmult(const std::vector< std::vector<int> > &, const int , const long int , const double* const , std::vector<double> &);
std::vector<double> dualtotals(const double *const, const int, const long int, const int *, const int );
int adjusttotals(std::vector<double> &, const double * const, const int, const long int, const int * const);
int makegreedy(std::vector<SuppPt> &, std::vector<int> & , int *, const int &, std::vector<double> &, const long int &);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   std::vector< std::time_t > pricingtimes;

   std::cout.precision(10);
   int Pnum =14;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   std::string startmonth = "Sep";
   int startyear = 14;

   int NumMoMeas = 1;//Number of Months in a single Measure. Use 1 for original approach of each month is its own measure
   //Otherwise, Number of measures that will be created. Will combine NumMoMeas of months into each measure
   int colmeth = 1;// 0: All potential columns, positive int: int number of best columns
   int probdel = 0;
   srand(time(NULL));
   
   int warmstart = 1; //If 1, the previous solution is saved and reintroduced before ->optimize()
   // This should be 1 for fastest running times
   
   int Psnum[Pnum];
//   int Psnumtemp[Pnum];
//   int numinprice = 2; // Numbers not = 1,2 have not been fully implemented
   std::vector<SuppPt> Psupp;
   int totalsupp = 0;
   int startindices[Pnum];
   int indices[Pnum];
   int endindices[Pnum];
//   int Aprn=0;
   double lambda[Pnum];
//   int warmstart = 1; //If 1, the previous solution is saved and reintroduced before ->optimize()
//   const double pricetol = -1e-6; // -1e-7 will usually allow duplication before termination. -1e-6 cuts off just a couple iterations with no apparent loss of accuracy in master objective

//   int comptoTrue = 0;
//   double Trueobj = 0;
//   std::ofstream outerror;
//   double currenttol = 1.0;
//
   std::string temp;

   std::ifstream indata;
   std::ostringstream fname;
   int year;
   std::string month;
   if (argc > 1)
   {
      fname << argv[1];
      indata.open(fname.str());
   }
   else
   {
      std::cout << "Please enter filepath for DenverCrime file." << std::endl;
      std::string filename;
      getline( std::cin, filename);
      indata.open(filename);
   }
   if (!indata.good())
   {
      std::cout << "File not found." << std::endl;
      return 1;
   }
   getline(indata, temp, '\n');// Advance past headers
   indata >> year;
   if (startyear < year)
   {
      std::cout << "Invalid Year Selection." << std::endl;
      return 1;
   }
   while (year < startyear)
   {
      getline(indata, temp, '\n');// Advance to next entry
      indata >> year;
   }
   indata >> month;
   while (month != startmonth)
   {
      getline(indata, temp, '\n');
      indata >> year;
      indata >> month;
      if (year != startyear)
      {
         std::cout << "Selected Month/Year not in data." << std::endl;
         return 1;
      }
      if (indata.eof())
      {
         indata.close();
         std::cout << "Invalid Input. End of file reached." << std::endl;
         return 1;
      }
   }
   
   double loc1;
   double loc2;
   indata >> loc1 >> loc2;
   std::string currentmonth = startmonth;
         
   for (int i = 0; i < Pnum; ++i)
   {
      Psnum[i] = 0;
      startindices[i] = totalsupp;
      indices[i] = totalsupp;
      for (int j = 0; j < NumMoMeas; ++j)
      {
         while (month == currentmonth )
         {
            //Adding this shift so all coordinates match the Pricing IP coordinates
            Psupp.push_back(SuppPt(loc1+115, loc2 , 1.0));
            ++Psnum[i];
            ++totalsupp;
                  
            indata >> year >> month >> loc1 >> loc2;
            if (indata.eof())
            {
               indata.close();
               std::cout << "Invalid Input. End of file reached while reading in months to measures." << std::endl;
               return 1;
            }
         }
         currentmonth = month;
      }
      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
      
      //Scale the masses for each month to sum to 1
      double totalmass = 0;
      for (int j = 0; j < Psnum[i]-1; ++j)
      {
         Psupp[startindices[i]+j].mass /= Psnum[i];
         totalmass += Psupp[startindices[i]+j].mass;
      }
      Psupp[totalsupp-1].mass = 1-totalmass;
      currentmonth = month;
      endindices[i] = totalsupp-1;
   }
   Psupp.resize(totalsupp);
   indata.close();
   
   //Compute the weights for each month
   double sum = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      lambda[i] = (double)Psnum[i]/totalsupp;
      sum += lambda[i];
   }
   lambda[Pnum-1] = 1-sum;

   
   int index = 0;
   // ensure P_i sum to exactly 1
   std::vector<SuppPt>::iterator Psupptest = Psupp.begin();
   for (int i = 0; i < Pnum; ++i)
   {
      validateDist(Psupptest, Psupptest+Psnum[i]);
      Psupptest += Psnum[i];
   }
   
   std::vector<SuppPt>::iterator Psuppit = Psupp.begin();
   for (int j = 0; j < totalsupp; ++j)
   {
      Psuppit->combos = j;
      //      std::cout << *Psuppit << " " << *Psupptempit << std::endl;
      ++Psuppit;
   }
   Psuppit = Psupp.begin();

   //Calculate number of possible supp pts S0
   long int S0 = 1;
   for (int i = 0; i < Pnum; ++i)
   {
      S0 *= Psnum[i];
   }
   std::cout << "Size of S0: " << S0 << std::endl;
   
   //For ensuring a column does not get added repeatedly
   std::vector<bool> colin(S0,false);
   //Temporary storage of the possible support points
   //In previous versions, all support points were stored
   SuppPt Pbar0;
   
   //Storing the vector of costs is necessary in this version of CG
   std::vector<double> c(S0,0);
   
   //Model for Master Problem
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   model->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   std::vector< GRBVar > w;
   std::vector< long int> wloc;
   GRBLinExpr  exp[totalsupp];
   
   //Compute the c's
   index = 0;
   for (unsigned long int j = 0; j < S0; ++j)
   {
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         sum1 += lambda[i]*Psupp[index].loc1;
         sum2 += lambda[i]*Psupp[index].loc2;
      }
         
      Pbar0 = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         c[j] += lambda[i]*((Pbar0.loc1-Psupp[index].loc1)*(Pbar0.loc1-Psupp[index].loc1) +(Pbar0.loc2-Psupp[index].loc2)*(Pbar0.loc2-Psupp[index].loc2));
      }
         
      int k = Pnum-1;
      if (indices[k] < endindices[k])
      {
         ++indices[k];
      }
      else
      {
         int temp = k-1;
         while (temp >= 0)
         {
            if (indices[temp] == endindices[temp])
            {
               temp -= 1;
            }
            else
            {
               ++indices[temp];
               break;
            }
         }
         for (int l = k; l > temp; --l)
         {
            indices[l] = startindices[l];
         }
      }
   }
   
   //Reset indices to start of each measure
   for (int i = 0; i < Pnum; ++i)
   {
      indices[i] = startindices[i];
   }
   
   //Create a temporary copy of the set of support points, will be modified
   std::vector<SuppPt> Psupp2(totalsupp);
   std::copy(Psuppit, Psuppit+totalsupp, Psupp2.begin() );
   std::vector<double> winit(S0,0);
   while (Psupp2[Psnum[0]-1].mass >1e-16)
   {
      //Find smallest remaining mass among support points of current combination
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[indices[i]].mass < minmass)
         {
            minmass = Psupp2[indices[i]].mass;
         }
      }
      
      //Index math based on Algorithm 1 in A Column Generation Approach to the Discrete Barycenter Problem
      long int denom = S0;
      int h = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[indices[i]].combos-index;
         denom /= Psnum[i];
         h += denom*index2;
         index += Psnum[i];
         Psupp2[indices[i]].mass -= minmass;
         if (Psupp2[indices[i]].mass <= 1e-16)
         {
            ++indices[i];
         }
      }
      
      winit[h] = minmass;
//      std::cout << h << std::endl;
   }
   
   long int k = 0;
   int jindex = 0;
   long int Pprod;
   
   std::vector<int> vbases;
   double oldobj;
   double newobj;
    
   //Introduce greedy-generated solution to initial master problem
   for (std::vector<double>::iterator wit = winit.begin(); wit != winit.end(); ++wit)
   {
         if (*wit != 0)
         {
            w.push_back(model->addVar(0.0, GRB_INFINITY, c[k], GRB_CONTINUOUS));//This updates the objective
            wloc.push_back(k);
            vbases.push_back(0);
            colin[k] = true;
            Pprod = S0/Psnum[0];
            jindex = floor(k/Pprod);
            exp[jindex] += *(w.end()-1);
            for (int l = 1; l < Pnum-1; ++l)
            {
               Pprod /= Psnum[l];
               jindex = startindices[l]+floor( (k % (Pprod*Psnum[l]))/Pprod);
               exp[jindex] += *(w.end()-1);
            }
            jindex = startindices[Pnum-1]+ k % Psnum[Pnum-1];
            exp[jindex] += *(w.end()-1);
         }
      ++k;
   }
   for (int j = 0; j < totalsupp; ++j)
   {
      model->addConstr(exp[j] == Psupp[j].mass);
   }
   
   t = std::chrono::steady_clock::now();
   std::cout << "Setup Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   std::cout << "Initial Number of Columns: " << w.size() <<std::endl;
   int iter = 1;
   
   model->optimize(); //First solve for dual values
   newobj = model->get(GRB_DoubleAttr_ObjVal);
   oldobj = newobj;
   auto t0 = t;
   t = std::chrono::steady_clock::now();
   std::cout << "First LP Solve Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-t0).count()  <<"ms" << std::endl;
 
   //Save solution information for warm start
   int cbases[totalsupp];
   GRBConstr* constrs =model->getConstrs();
   if (warmstart == 1)
   {
      for (int i = 0; i < w.size(); ++i)
      {
         vbases[i] = (w[i].get(GRB_IntAttr_VBasis));
      }
      for (int i = 0; i < totalsupp; ++i)
      {
         cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
      }
   }
   
   std::vector<int> vbstart = vbases;
   bool isneg = true;
   long int maxcolnum = w.size();
   int maxiter = 25;
   int total_iter = 0;
   while ( isneg )
   {
      if (total_iter >= maxiter)
      {
         break;
      }
      ++total_iter;
      isneg = false;
//      double currentobj = model->get(GRB_DoubleAttr_ObjVal);
//      std::cout << currentobj << std::endl;
      t0 = t;
      t = std::chrono::steady_clock::now();
      constrs =model->getConstrs();
      double* yhat = model->get(GRB_DoubleAttr_Pi, constrs, totalsupp);

      std::vector< GRBVar >::iterator wit = w.begin();
      long int vloc = 0;
      if (newobj != oldobj)
      {
         if (w.size() > maxcolnum)
         {
            maxcolnum = w.size();
         }

         if (probdel > 0)
         {
            while( wit < w.end())
            {
//               std::cout << "Index is: " << vloc;
//               std::cout << " vbases value is: " << vbases[vloc];
//               std::cout << " Number of variables is: " <<w.size() <<std::endl;
               if (vbases[vloc] == -1 && (rand() %100) < probdel)
               {
                  //remove variable from LP
                  model->remove(w[vloc]);
                  //delete variable from list of variables w
                  w.erase(wit);
                  //set toggle to not included
      //                  std::cout << wloc[vloc] <<std::endl;
                  colin[wloc[vloc]] = false;
                  wloc.erase(wloc.begin()+vloc);
                  //delete from vbases for warm start

                  vbases.erase(vbases.begin()+vloc);
               }
               else
               {
                  ++vloc;
                  ++wit;
               }
            }
            model->update();
            vbstart = vbases;
         }
      }
//      std::cout << "Number of variables is " << w.size() << " and in vbasis " << vbases.size() <<std::endl;
      std::vector<double> yhatA= dualtotals( yhat,  totalsupp, S0, Psnum, Pnum);

      std::vector< double> min(colmeth,0);
      std::vector<long int> minloc(colmeth,0);

      //This method is the naive, just add all of them that are sufficiently negative back in.
      if (colmeth == 0)
      {
         for (int j = 0; j < S0; ++j)
         {
            if (c[j] - yhatA[j]< -1e-14 && colin[j] == false)
            {
               isneg = true;
               double expcoeff2[totalsupp];
               for (int l = 0; l < totalsupp; ++l)
               {
                  expcoeff2[l] = 0;
               }
               Pprod = S0/Psnum[0];
               jindex = floor(j/Pprod);
               expcoeff2[jindex] = 1;
               for (int l = 1; l < Pnum-1; ++l)
               {
                  Pprod /= Psnum[l];
                  jindex = startindices[l]+floor( (j % (Pprod*Psnum[l]))/Pprod);
                  expcoeff2[jindex] =1;
               }
               jindex = startindices[Pnum-1]+ j % Psnum[Pnum-1];
               expcoeff2[jindex] = 1;
               GRBColumn newcol;
               newcol.addTerms(expcoeff2, constrs, totalsupp);
               w.push_back(model->addVar(0.0, GRB_INFINITY, c[j], GRB_CONTINUOUS, newcol));//This updates the objective
               wloc.push_back(j);
               vbases.push_back(-1);
               colin[j] = true;
            }
         }
      }
      else
      {
         //will introduce colmeth number of columns
         for (int j = 0; j < S0; ++j)
         {
            double cnew = c[j]-yhatA[j];
            if (cnew < min.back() && colin[j] == false)
            {
               std::vector<double>::iterator minit = min.end()-1;
               std::vector<long int>::iterator minlocit = minloc.end()-1;
//               std::cout << "Start of search with value: " << cnew << std::endl;
               for (int i = colmeth-1; i >= 0; --i)
               {
                  if (cnew > *minit )
                  {
                     min.insert(minit+1,cnew);
                     isneg = true;
                     minloc.insert(minlocit+1,j);
                     min.pop_back();
                     minloc.pop_back();
                     break;
                  }
                  else if (minit == min.begin())
                  {
                     min.insert(minit,cnew);
                     isneg = true;
                     minloc.insert(minlocit,j);
                     min.pop_back();
                     minloc.pop_back();
                  }
                  else
                  {
                     --minit;
                     --minlocit;
                  }
               }
            }
         }
         for (std::vector<double>::iterator minit = min.end()-1; minit >= min.begin(); --minit)
         {
            if (*minit == 0)
            {
               min.pop_back();
               minloc.pop_back();
            }
         }
//         std::cout << "After removing extra 0's, min is: ";
//         for (std::vector<double>::iterator minit = min.begin(); minit != min.end(); ++minit)
//         {
//            std::cout << *minit << " ";
//         }
//         std::cout << std::endl;
      }
      t0 = t;
      t = std::chrono::steady_clock::now();
      pricingtimes.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t-t0).count() );
      if (isneg)
      {
         if (colmeth >0)
         {
            for (int i = 0; i < min.size(); ++i)
            {
               double expcoeff2[totalsupp];
               for (int l = 0; l < totalsupp; ++l)
               {
                  expcoeff2[l] = 0;
               }
               Pprod = S0/Psnum[0];
               jindex = floor(minloc[i]/Pprod);
               expcoeff2[jindex] = 1;
               for (int l = 1; l < Pnum-1; ++l)
               {
                  Pprod /= Psnum[l];
                  jindex = startindices[l]+floor( (minloc[i] % (Pprod*Psnum[l]))/Pprod);
                  expcoeff2[jindex] =1;
               }
               jindex = startindices[Pnum-1]+ minloc[i] % Psnum[Pnum-1];
               expcoeff2[jindex] = 1;
               GRBColumn newcol;
               newcol.addTerms(expcoeff2, constrs, totalsupp);
               w.push_back(model->addVar(0.0, GRB_INFINITY, c[minloc[i]], GRB_CONTINUOUS, newcol));//This updates the objective
               wloc.push_back(minloc[i]);
               vbases.push_back(-1);
               colin[minloc[i]] = true;
            }
         }
         if (warmstart == 1 )
         {
            model->update();
            for (int i = 0; i < w.size(); ++i)
            {
               w[i].set(GRB_IntAttr_VBasis, vbases[i]);
            }
            for (int i = 0; i < totalsupp; ++i)
            {
               constrs[i].set(GRB_IntAttr_CBasis, cbases[i]);
            }
         }
         vbstart = vbases;
         oldobj = newobj;
         model->optimize();
//         std::ostringstream fname2;
//         fname2 <<"/Users/spatterson/ForXcode/Barycenter/Column Generation/TestNewGen_"<<iter<<".lp";
//         model->write(fname2.str());
         newobj = model->get(GRB_DoubleAttr_ObjVal);
         ++iter;
         if (warmstart == 1)
         {
            for (int i = 0; i < w.size(); ++i)
            {
               vbases[i] = (w[i].get(GRB_IntAttr_VBasis));
//               std::cout << vbases[i] <<std::endl;
            }
            for (int i = 0; i < totalsupp; ++i)
            {
               cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
            }
         }
      }
   }
   
   t = std::chrono::steady_clock::now();
   std::cout << "Total Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   if (probdel > 0)
   {
      std::cout << "Maximum Column Number: " << maxcolnum << std::endl;
   }
   std::cout << "Total Columns Introduced: " << w.size() << std::endl;
   std::cout << "Optimal value: " << newobj << " found after " << iter << " solves." << std::endl;
   std:: cout << "Average Pricing Time " << std::accumulate(pricingtimes.begin(), pricingtimes.end(),0)/iter << "ms" <<std::endl;
   
   std::cout << "Enter to end." <<std::endl;
   std::cin.get();
   model->terminate();
   model->reset();
   delete model;
   delete env;
   return 0;
}

   bool validateDist(const std::vector<SuppPt>::iterator it1, const std::vector<SuppPt>::iterator it2)
   {
      std::vector<SuppPt>::iterator it = it1;
      double total = 0;
      while (it != it2)
      {
         total += it->mass;
         ++it;
      }
      
      if (total < 0.95 || total > 1.05)
      {
         std::cout << "Warning: Far from 1. Consider alternative." << total << std::endl;
      }
      if (total != 1)
      {
         it1->mass = it1->mass + (1.0-total);
      }
      return true;

}

int makegreedy(std::vector<SuppPt> &Psupp2, std::vector<int> &Psnum2, int * Psnum, const int &Pnum, std::vector<double> &winit,const long int &S0)
{
   std::vector<int> indices(Pnum,0);
   std::vector<int> startindices(Pnum,0);
   
   for (int i = 1; i < Pnum; ++i)
   {
      indices[i] = indices[i-1];
      indices[i] += Psnum2[i-1];
      startindices[i] = indices[i];
   }
   
   while (Psupp2[Psnum2[0]-1].mass >1e-12)
   {
      //      std::cout << Psupp2[Psnum2[0]-1].mass <<std::endl;
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[indices[i]].mass < minmass)
         {
            minmass = Psupp2[indices[i]].mass;
         }
      }
      
      //Index math based on Algorithm 1 in A Column Generation Approach to the Discrete Barycenter Problem
      long int denom = S0;
      int h = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[indices[i]].combos-index;
         denom /= Psnum[i];
         h += denom*index2;
         index += Psnum[i];
         Psupp2[indices[i]].mass -= minmass;
         if (Psupp2[indices[i]].mass <= 1e-12)
         {
            ++indices[i];
         }
      }
      
      winit[h] = minmass;
   }
   return 0;
}


std::vector<double> dualtotals(const double * const yhat, const int Amrn, const long int blocklength, const int * Psnum, const int Pnum)
{
   std::vector<double> ytotals(blocklength,0);
   int i =0; //Should really pass this in. Related to numinprice
   int totalPs = 0;
   long int startindex = 0;
   long int blockinum = 1;
   long int blockilength = blocklength;
   long int blocki1length = blockilength/Psnum[i];
   long int j = 0;
   long int index = 0;
   while (j < Amrn)
   {
      if (yhat[j] != 0)
      {
         //         std::cout << yhat[j] << std::endl;
         index = 0;
         for (int k = 0; k < blockinum; ++k)
         {
            for (int l = 0; l < blocki1length; ++l)
            {
               //               std::cout << "index: " << startindex+index << " yhat: " << yhat[j] << std::endl;
               ytotals[startindex+index] += yhat[j];
               ++index;
            }
            index = blockilength*(k+1);
         }
      }
      ++j;
      //      std::cout << j << std::endl;
      startindex += blocki1length;
      if (j-totalPs == Psnum[i])
      {
         blockinum *= Psnum[i];
         blockilength /= Psnum[i];
         startindex = 0;
         totalPs += Psnum[i];
         ++i;
         if (i < Pnum)
         {
            blocki1length /= Psnum[i];
         }
      }
   }
   return ytotals;
}
