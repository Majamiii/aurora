#define WIN32_LEAN_AND_MEAN 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <cmath>
#include <time.h>
#include <thread>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */

//my own stuff
#include "utility_functions.h"
#include "simulation.h"
#include "electrostatics.cpp"
#include "io.cpp"
#include "electron.h"

#include "utility_functions.h"

#include <windows.h>

#include <ctime>
#include <random>
#include <cstdlib>

#include <omp.h>
// OMP_NUM_THREADS
#ifndef THREAD_NUM
#define THREAD_NUM 12
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;


int main() {
   omp_set_num_threads(THREAD_NUM);
   srand ( time(NULL));
   vector<float> randoms;
   Simulation sim;
   sim.init();
   PhotonDensity photon_density;
   photon_density.init(&sim);
   EnergyDensity energy_density;
   energy_density.init(photon_density.resolution_z);
   MagneticField mag_field;
   ChargeDensity rho;
   ElectricField E_field;
   E_field.init(&sim, &rho);
   
   int voxelx,voxely,voxelz;
   int Evoxelx,Evoxely,Evoxelz;
   vector<float> B;
   float By=-1.0;
   float Bx=0.0;
   float Bz=0.0;
   float tmpfloat;

   double verovatnoca_arr[359];
   int ver_i = 0;
   double z_2nd = 0;

   float v_curr =  0;
   float alpha_0 = 0;

   float suma_verovatnoca = 0;

   srand (static_cast <unsigned> (time(0)));


   for (int i=0; i<sim.N; i++) {
      Electron ee;
      ee.sim = &sim;
      ee.reset();
      ee.ID = i;
      sim.electrons.push_back(ee);

   }
   cout << "Electron objects created" << endl;
 
   ofstream debug_xyz;

   float rnd;
   float p_interaction = 0.1;
   int t=0;
   // int t0=0;
   cout << endl;
   cout << "Beginning to integrate: " << endl;
   for (t=0; t< sim.tmax; t++ ) {
      cout<<" t : "<< t << "\n";

     sim.t = t;
      if (t%5==0 && t>19) {
         photon_density.write_image(t);
         photon_density.reset();
         energy_density.write_out(t);
         energy_density.reset();
      }

      // #pragma omp parallel for private(voxelx, voxely, voxelz, suma_verovatnoca) 

      for (int i=0; i<sim.N ; i++) {      //   for every electron

         Electron* e = &sim.electrons[i]; //make a pointer to the electron we're dealing with
         
         //if electron has left the cell, or has not enough energy remaining, we're done tracking it.
         //Reset it.
         while (  
            e->z!=e->z //(that will evaluate true if e->z == nan) 
            || e->z < 0 
            || e->z >= 1.05*(sim.box_sizez)
            || e->E <= (sim.hc / e->sim->wavelength_red )  
          ) {
            e->reset();
            e->t = t;
         }
          
         voxelx = int((float)photon_density.resolution_x / (float)sim.box_sizex * e->x);
         voxely = int((float)photon_density.resolution_y / (float)sim.box_sizey * e->y);
         voxelz = int((float)photon_density.resolution_z / (float)sim.box_sizez * e->z);
         Evoxelx = voxelx;
         Evoxely = voxely;
         Evoxelz = voxelz;
        
         e->calculate_probabilities(e->z);

         rnd = (float)rand()/RAND_MAX;
         if (rnd/3.0 < e->p_emit) {
            e->interaction_count++;
            rnd = (float)rand()/RAND_MAX;
            if (rnd < e->p_emit_r) {      //kiseonik
               //emit red
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_red();
               e->emitting_wavelength = sim.wavelength_red;
               z_2nd = 8*8;      //atomic num of oxygen, for the rutherford formula
            } else if (rnd < (e->p_emit_r + e->p_emit_g)) {    //kiseonik
               //emit green
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_green();
               e->emitting_wavelength = sim.wavelength_green;
               z_2nd = 8*8;      //atomic num of oxygen, for the rutherford formula

            } else {                      //azot
               //emit blue
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_blue();
               e->emitting_wavelength = sim.wavelength_blue;
               z_2nd = 7*7;      //atomic num of nitrogen, for the rutherford formula
            }

         }

         //Check to see if electron will emit energy this timestep:
         if (e->emitting==1 && e->dead_counter == 0) {
            e->E -= sim.E_loss_factor*sim.hc / e->emitting_wavelength;
            energy_density.increment(voxelz,sim.E_loss_factor*sim.hc / e->emitting_wavelength); //deltaE = hc/lambda
            e->emitting_time_left -= 1;
            
            if (  e->z < sim.box_sizez
                  && e->z >=0
                  && e->x >= 0 
                  && e->x <= sim.box_sizex
                  && e->y >= 0
                  && e->y <= sim.box_sizey    ) {

               if (e->emitting_wavelength == sim.wavelength_red ) {        //deltaE = 2479,68
                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.8) ;
                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.16) ;
                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.4) ;
               } else if (e->emitting_wavelength == sim.wavelength_green ) {//deltaE = 2066
                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.18) ;
                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.69) ;
                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.4) ;
               } else if (e->emitting_wavelength == sim.wavelength_blue ) {//deltaE=1549
                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.4) ;
                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.4) ;
                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.2) ;
               }

            }

            if (e->emitting_time_left < 1) { //if it's done emitting
               e->emitting = 0;
               e->emitting_wavelength = 0;
            }
         }

         


//CHARGE DENSITY CALCULATION///////////////////////////////
        voxelx = int((float)rho.resolution_x / (float)sim.box_sizex * e->x); 
        voxely = int((float)rho.resolution_y / (float)sim.box_sizey * e->y);
        voxelz = int((float)rho.resolution_z / (float)sim.box_sizez * e->z);
        if (
             voxelx >= 0 
          && voxelx<rho.resolution_x 
          && voxely >= 0 
          && voxely<rho.resolution_y 
          && voxelz >= 0 
          && voxelz<rho.resolution_z) 
        {
            rho.incr_element(voxelx,voxely,voxelz);  
        }
///////////////////////////////////////

         
//FORCES////////////////////////
         e->Fx = 0;
         e->Fy = 0;
         e->Fz = 0;
         //magnetic force F = v x B
         e->Fx += (e->vy * Bz - e->vz * By) / 1e-1 ;
         e->Fy += (e->vz * Bx - e->vx * Bz) / 1e-1;
         e->Fz += (e->vx * By - e->vy * Bx) / 1e-1 ;
         //electric force F = q E
         e->Fx += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ex,voxelx,voxely,voxelz));
         e->Fy += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ey,voxelx,voxely,voxelz));
         e->Fz += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ez,voxelx,voxely,voxelz));
////////////////////////////////




   // SCATTERINGGGGGGGGGGG

         int angle=0;
         float k=0;

         int scale = 0;

         float sum_cdf = 0;

         float radian=0;
         float sin_4th=0;
         float cnst = 3.4e-5; 

         float verovatnoca = 0;

         for(angle=0; angle<180; angle++) {

            radian=(angle)*M_PI/180;
            sin_4th = pow(sin(radian), 4);

            if (sin(radian) == 0) {
               sin_4th = pow((sin((angle+1)*M_PI/180)),4);
            }

            // if(angle>180) {
            //    sin_4th = -sin_4th;
            // }

            k = cnst / sin_4th;
            verovatnoca = k / (pow(e->E,2)) * z_2nd; 

            suma_verovatnoca += verovatnoca;       //za CDF

            verovatnoca_arr[angle] = verovatnoca;

         }

         for(scale=1; scale<180; scale++) {
            // suma_verovatnoca += verovatnoca_arr[scale];
            verovatnoca_arr[scale] += verovatnoca_arr[scale-1];
         }

         for(scale=0; scale<180; scale++) {
            verovatnoca_arr[scale] = verovatnoca_arr[scale] / suma_verovatnoca;
            // cout<<"\n verovatnoca: "<< verovatnoca_arr[scale] << "   za ugao  "<<scale;
         }
         suma_verovatnoca = 0;
         sum_cdf = 0;



//INTEGRATION////////////////////////////////////
         // Get a different random number each time the program runs

         float randomNum = 0 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(1)));
         float min_dif = 100;
         float dif = 0;
         float alpha = 0;
         int find_c = 0;

         for(find_c = 0; find_c < 180; find_c++) {
            dif = abs(verovatnoca_arr[find_c]-randomNum);
            if (dif<min_dif) {
               min_dif = dif;
               alpha = find_c;
            }
         }

         int randomNum2 = rand() % 2;

         if ((randomNum2 == 1)&&(alpha != 0)) {
            alpha = -alpha;
         }

         // cout << "\n angle: "<<alpha;

         min_dif = 100;

         if((abs(e->vx)>0.01)&&(abs(e->vy)>0.01)) {            
            alpha_0 = atan(e->vy / e->vx);
            alpha += alpha_0;
         }
         else {
            alpha_0 = 0;
         }


         // cout<<"\n angle0: "<<alpha;
         v_curr = sqrt(pow(e->vx, 2) + pow(e->vy, 2));

         //Perform equation of motion integration:
         e->vx +=  (e->Fx / sim.m_e) * sim.dt ;

         e->vx = v_curr * cos(alpha*3.14159/180);         //converting to radians

         e->vy +=  (e->Fy / sim.m_e * sim.dt );
         e->vy = v_curr * sin(alpha*3.14159/180);

         e->vz +=  (e->Fz / sim.m_e) * sim.dt ;

         e->x += e->vx * sim.dt;
         e->y += e->vy * sim.dt;
         e->z += e->vz * sim.dt;
//////////////////////////////////////////////////
         

///////////PERIODIC BOUNDARY CONDITIONS/////
           while (e->x < 0) {e->x += sim.box_sizex;}
           while (e->y < 0) {e->y += sim.box_sizey;}
           while (e->x > sim.box_sizex) {e->x -= sim.box_sizex;}
           while (e->y > sim.box_sizey) {e->y -= sim.box_sizey;}
///////////////////////////////////////////


   // RESCALING VELOCITIES

         // float R = 0.5*sim.m_e*(e->vx*e->vx + e->vy*e->vy + e->vz*e->vz) / e->E;
         // e->vz = e->vz/sqrt(R);

   // RESCALING VELOCITIES



// //print out a bunch of info for one of the particles so we can see how simulation is progressing 
         if (e->ID==3) { 
            cout << t 
                 << "\t" 
                 << e->x 
                 << "\t" 
                 << e->y 
                 << "\t" 
                 << e->z 
                 << "\tE: " 
                 << e->E 
                 <<  "\t"  
                 << e->Fx 
                 << "\t" 
                 << e->Fy 
                 << "\t" 
                 << e->Fz 
                 << "\t" 
                 << e->vx 
                 << "\t" 
                 << "  y: "<<e->vy 
                 << "\t" 
                 <<"   z: "<<e->vz 
                 << "\n" ;
         }
         

      } // end of loop over electrons

      //recompute the electric field
      if (t%sim.E_field_recalc==0 ) {
         E_field.compute();
      }
      rho.reset(); //zero out the charge density for next time

   } //end of loop over time
   
   cout << endl;

   cout << "\nCLEANING UP" << endl;

   energy_density.cleanup();
   photon_density.cleanup();
   E_field.cleanup();
   rho.cleanup();

   cout << "\n\n" << endl;

   return 0;

}








