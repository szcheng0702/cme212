#include <iostream>
#include <cassert>
#include <ctime>
#include <unistd.h>
#include <vector>

#include <boost/multi_array.hpp>

#include "Viewer.hpp"

 /** 
 * @class bz_model
 * @brief A class to demonstrate Belousov-Zhabotinsky Reaction.
 *
 * Users can define the grid size, number of reactants, 
 * and the box size used in the diffusion model */
template<typename I, typename T,typename array_type>
class bz_model {
	using viewer_type = Viewer<I, T>;

public:
	bz_model(I xsize, I ysize, I num_reactants, I boxsizem, I boxsizen, array_type& y1)
		: xsize_(xsize),
		ysize_(ysize),
		num_reactants_(num_reactants),
		m((boxsizem-1)/2),
		n((boxsizen-1)/2),
		y(y1)
	{}

	I getXsize()
	{
		return xsize_;
	}

	I getYsize()
	{
		return ysize_;
	}

	/** Helper function to perform diffusion process
	*  @pre @a shape(v)=={@a xsize_,@a ysize_, @a num_reactants_}
	*  @pre for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, @ 0 <k< @a num_reactants_, @a y[i][j][k]>0
	*  @pre for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, sum(@a y[i][j][k])=1 if sum over  @ 0 <k< @a num_reactants_
	*  
	*  @post for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, sum(@a y_new[i][j][k])=1 if sum over  @ 0 <k< @a num_reactants_
	*  @post for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, @ 0 <k< @a num_reactants_, @a y_new[i][j][k]>0
	*/
	void updatediffusion(){
	//diffusion species update
		for (I i = 0; i < xsize_; ++i) {
			for (I j = 0; j < ysize_; ++j) {
				
				std::vector<T> summation (num_reactants_,0);
				for (I p = i-m; p <= i+m; ++p) {
					for (I q = j-n; q <= j+n; ++q) {
						//check if index exceeds the left boundary and get an updated index
						I temp1 = p < 0 ? (xsize_ + p) : p;
						//check if index exceeds the right boundary and get an updated index
						I id1 = temp1 > (xsize_ - 1) ? (temp1-xsize_) : temp1;

						//check if index exceeds the upper boundary and get an updated index
						I temp2 = q < 0 ? (ysize_ + q) : q;
						//check if index exceeds the lower boundary and get an updated index
						I id2 = temp2 >(ysize_ - 1) ? (temp2 - ysize_) : temp2;
		
						for(I k=0;k<num_reactants_;++k){
							summation[k]+=y[id1][id2][k];
						}

					}
				}

				//normalization
				for(I k=0;k<num_reactants_;++k){
					y[i][j][k] = 1.0 / ((2 * m + 1)*(2 * n + 1))*summation[k];
				}

			}
		}

	}

	/** Helper function to perform reaction process
	*  @pre @a shape(v)=={@a xsize_,@a ysize_, @a num_reactants_}
	*  @pre for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, @ 0 <k< @a num_reactants_, @a y[i][j][k]>0
	*  @pre for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, sum(@a y[i][j][k])=1 if sum over  @ 0 <k< @a num_reactants_
	*  
	*  @post for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, sum(@a y_new[i][j][k])=1 if sum over  @ 0 <k< @a num_reactants_
	*  @post for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, @ 0 <k< @a num_reactants_, @a y_new[i][j][k]>0
	*/
	void updatereaction() {
		//reaction in cell update
		for (I i = 0; i < xsize_; ++i) {
			for (I j = 0; j < ysize_; ++j) {
				//Has to hard code here, since I have no idea what the general equation is like
				//when number of reactants is larger than 3

				//update the mass fraction during reaction process, and set the floor to be 0
				T ya_temp=std::max(0.0,y[i][j][0]+y[i][j][0]*(ka*y[i][j][1]-kc*y[i][j][2]));
				T yb_temp=std::max(0.0,y[i][j][1]+y[i][j][1]*(kb*y[i][j][2]-ka*y[i][j][0]));
				T yc_temp=std::max(0.0,y[i][j][2]+y[i][j][2]*(kc*y[i][j][0]-kb*y[i][j][1]));

				//normalization
				T temp_sum = ya_temp + yb_temp + yc_temp;
				y[i][j][0] = ya_temp / temp_sum;
				y[i][j][1] = yb_temp/ temp_sum;
				y[i][j][2] = yc_temp / temp_sum;
			}
		}


	}


	/** Helper function to update the SFML viewer and color the cells
	*  @pre @a shape(v)=={@a xsize_,@a ysize_, @a num_reactants_}
	*  @pre for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, @ 0 <k< @a num_reactants_, @a y[i][j][k]>0
	*  @pre for any 0 <= i < @a xsize_, @ 0 < j < @a ysize_, sum(@a y[i][j][k])=1 if sum over  @ 0 <k< @a num_reactants_
	*  @pre @a num_reactants_=3
	*/
	void colorCells(viewer_type& v)
	{
		for (I i = 0; i < xsize_; ++i)
			for (I j = 0; j < ysize_; ++j)
				//Here the number of reactants has to be 3 to make the v.updatecell work
				v.updateCell(i, j, y[i][j][0], y[i][j][1], y[i][j][2]);
		v.refresh();
	}

	
private:
	I xsize_;
	I ysize_;
	I num_reactants_;
	I m;//boxsize in x direction
	I n;//boxsize in y direction
	array_type y;
	double ka=1;
	double kb=1;
	double kc=1;
};

int main()
{
  using index_type = int;
  using real_type  = double;
  using array_type = boost::multi_array<real_type,3>;
  
  std::string windowTitle = "Belousov-Zhabotinsky Reaction";
  index_type xsize = 100;
  index_type ysize = 100;
  index_type num_reactants=3;


  
  //Initialization
  srand((unsigned int)time(0));

  array_type y(boost::extents[xsize][ysize][num_reactants]);
  
  for (index_type i = 0; i < xsize; ++i) {
	  for (index_type j = 0; j < ysize; ++j) {
	  	std::vector<real_type> initialize_value(num_reactants);

	  	//We keep track of the sum because we want to do the initialization
	  	real_type thesum=0;


	  	for (index_type k=0;k<num_reactants;++k) {
	  	//randomly generate ya[i][j] and yb[i][j] with value less than 0.5
		  initialize_value[k] = ((double)rand() / (RAND_MAX));
		  thesum+=initialize_value[k];
		}
		
		/** We do this normalization because we want to keep the initialized value
         * satisfying ya+yb+yc=1 */
		for (index_type k=0;k<num_reactants;++k) {
			y[i][j][k]=initialize_value[k]/thesum;

		}
	 }
  }


  // Create a bz_model object
  bz_model<index_type, real_type,array_type> bz(xsize, ysize,3,3,3,y);
  
  // Create xsize by ysize display
  Viewer<index_type, real_type> v(xsize, ysize, windowTitle);
  
  while (v.isOpen())
  {
    v.event_loop();

    // update 
	bz.updatediffusion();
	bz.updatereaction();

	//Color the update cells
    bz.colorCells(v);

   
	//usleep(1000000); 

  }
  return 0;
}

