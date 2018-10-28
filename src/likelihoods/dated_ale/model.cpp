#include "exODT.h"
using namespace std;
using namespace bpp;

#include <bitset>


//p(ale) calculates Pi(Gamma) cf. ALEPAPER
scalar_type exODT_model::p(approx_posterior *ale)
{
  ale_pointer=ale;

  for (std::map<long int, std::map< scalar_type, std::map<int, scalar_type> > >::iterator it=q.begin();it!=q.end();it++)
    {
      for ( std::map< scalar_type, std::map<int, scalar_type> >::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	(*jt).second.clear();
      (*it).second.clear();
    }
  q.clear();

  //directed partitions and their sizes
  vector <long int>  g_ids;//del-loc
  vector <long int>  g_id_sizes;//del-loc
  for (map <int, vector <long int > > :: iterator it = ale->size_ordered_bips.begin(); it != ale->size_ordered_bips.end(); it++)
    for (vector <long int >  :: iterator jt = (*it).second.begin(); jt != (*it).second.end(); jt++)
      {
	g_ids.push_back((*jt));
	g_id_sizes.push_back((*it).first);
      }
  //root bipartition needs to be handled separately
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);

    // gene<->species mapping
  for (int i=0;i<(int)g_ids.size();i++)
    {
      long int g_id=g_ids[i];
      for (int rank=0;rank<last_rank;rank++)
	{
	  int n=time_slices[rank].size();
	  for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	    {
	      scalar_type t=time_slice_times[rank][t_i];
	      for (int branch_i=0;branch_i<n;branch_i++)
		{
		  int e = time_slices[rank][branch_i];
		  q[g_id][t][e]=0;
		}
	      q[g_id][t][alpha]=0;
	    }
	}

      if (g_id_sizes[i]==1)
	{
     /*   int id = 0;
        boost::dynamic_bitset<> temp = ale->id_sets[g_id];
        for (auto i = 0; i < ale->Gamma_size + 1; ++i) {
           // if ( BipartitionTools::testBit ( temp, i) ) {
                if ( temp[i] ) {
                    id = i;
                break;
            }
        }*/
        int id = 0;
        for (auto i=0; i< ale->Gamma_size + 1; ++i) {
            if ( ale->id_sets[g_id][i] ) {
                id=i;
                break;
            }
        }

        string gene_name=ale->id_leaves[ id /*g_id*/ ];
	//  string gene_name=ale->id_leaves[ (* (ale->id_sets[g_id].begin()) )];
	  vector <string> tokens;
	  boost::split(tokens,gene_name,boost::is_any_of(string_parameter["gene_name_separators"]),boost::token_compress_on);
	  string species_name;
	  if ((int)scalar_parameter["species_field"]==-1)
	    species_name=tokens[tokens.size()-1];
	  else
	    species_name=tokens[(int)scalar_parameter["species_field"]];
	  gid_sps[g_id]=species_name;
	}
    }

  for (int i=0;i<(int)g_ids.size();i++)
    {
      // directed partition (dip) gamma's id
      bool is_a_leaf=false;
      long int g_id=g_ids[i];
      if (g_id_sizes[i]==1)
	is_a_leaf=true;

      vector <long int> gp_ids;//del-loc
      vector <long int> gpp_ids;//del-loc
      vector <scalar_type> p_part;//del-loc
      if (g_id!=-1)
	for (unordered_map< pair<long int, long int>,scalar_type> :: iterator kt = ale->Dip_counts[g_id].begin(); kt != ale->Dip_counts[g_id].end(); kt++)
	  {
	    pair<long int, long int> parts = (*kt).first;
	    long int gp_id=parts.first;
	    long int gpp_id=parts.second;
	    gp_ids.push_back(gp_id);
	    gpp_ids.push_back(gpp_id);
	    if (ale->Bip_counts[g_id]<=scalar_parameter["min_bip_count"])
	      p_part.push_back(0);
	    else
	      p_part.push_back(ale->p_dip(g_id,gp_id,gpp_id));
	  }
      else
	{
	  //root bipartition needs to be handled separately
        map<set<long int>,int> bip_parts;
        for (map <long int,scalar_type> :: iterator it = ale->Bip_counts.begin(); it != ale->Bip_counts.end(); it++)
	    {
            long int gp_id=(*it).first;
            boost::dynamic_bitset<> gamma =ale->id_sets.at(gp_id);
            boost::dynamic_bitset<> not_gamma = ~gamma;
            not_gamma[0] = 0;
            long int gpp_id = ale->set_ids.at(not_gamma);

            set <long int> parts;
            parts.insert(gp_id);
            parts.insert(gpp_id);
            bip_parts[parts]=1;
            // gamma.clear();
            // not_gamma.clear();
	    }
	  for (map<set<long int>,int> :: iterator kt = bip_parts.begin();kt!=bip_parts.end();kt++)
	    {
	      vector <long int> parts;
            for (set<long int>::iterator sit=(*kt).first.begin();sit!=(*kt).first.end();sit++) {
                parts.push_back((*sit));
            }
	      long int gp_id=parts[0];
	      long int gpp_id=parts[1];
	      gp_ids.push_back(gp_id);
	      gpp_ids.push_back(gpp_id);

            //Here we can create a new ale->Bip_counts[gp_id], in particular for leaves.
            //We may want to add the leaf entries for Bip_counts when Bip_counts is first created.
	      if (ale->Bip_counts[gp_id]<=scalar_parameter.at("min_bip_count") and not ale->Gamma_size<4)
		p_part.push_back(0);
	      else
		p_part.push_back(ale->p_bip(gp_id));
	    }
	  bip_parts.clear();
	}
      int N_parts=gp_ids.size();

      //iterate over all postions along S
      for (int rank=0;rank<last_rank;rank++)
	{
	  int n=time_slices[rank].size();
	  for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	    {
	      //######################################################################################################################
	      //#########################################INNNER LOOP##################################################################
	      //######################################################################################################################

	      scalar_type t=time_slice_times[rank][t_i];
	      scalar_type tpdt,tpdt_nl;
	      if ( t_i < (int)time_slice_times[rank].size()-1 )
		tpdt=time_slice_times[rank][t_i+1];
	      else if (rank<last_rank-1)
		tpdt=time_slice_times[rank+1][0];
	      else
		//top of root stem
		tpdt=t_begin[time_slices[rank][0]];

	      if (scalar_parameter["event_node"]==1 and false)
		tpdt_nl=t;
	      else
		tpdt_nl=tpdt;

	      //root
	      scalar_type Delta_t=(tpdt-t)*1;

	      //Delat_bar corresponds to sigma in ALEPAPER
	      scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];
	      //scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank];
	      //scalar_type tmp
	      scalar_type p_Delta_bar=Delta_bar*Delta_t;
	      scalar_type Ebar=Ee[-1][t];

	      //boundaries for branch alpha virtual branch

	      //boundary at present
	      if (t==0)
		q[g_id][t][alpha]=0;

	      //boundary between slice rank and rank-1 slice is trivial
	      ;//q[g_id][t][alpha]=q[g_id][t][alpha];

	      //boundaries for branch alpha virtual branch.
	      if(1)
		{
		  for (int branch_i=0;branch_i<n;branch_i++)
		    {
		      int e = time_slices[rank][branch_i];

		      //boundaries for branch e
		      //boundary at present
		      if (t==0)
			{
                if (is_a_leaf && extant_species[e]==gid_sps[g_id])	{
			    q[g_id][t][e]=1;
                }
			  else
			    q[g_id][t][e]=0;
			}
		      //boundary between slice rank and rank-1
		      else if (t_i==0)
			{
			  //terminating branch is last in time_slices and defines a represented speciation
			  if (branch_i==n-1 && rank>0)
			    {
			      int f=daughters[e][0];
			      int g=daughters[e][1];
			      scalar_type Eft=Ee[f][t];
			      scalar_type Egt=Ee[g][t];

			      scalar_type q_sum=0;
			      //q[g_id][t][e]=0;

			      scalar_type SL_fLg=q[g_id][t][f]*Egt;
			      scalar_type SL_Lfg=q[g_id][t][g]*Eft;
			      //SL EVENT, events #3 and #4 in part c of Fig.A1 in http://arxiv.org/abs/1211.4606
			      //q[g_id][t][e]=q[g_id][t][f]*Egt + q[g_id][t][g]*Eft;
			      q_sum+=SL_fLg+SL_Lfg;
			      //SL.

			      //non-leaf directed partition
			      if (not is_a_leaf)
				for (int i=0;i<N_parts;i++)
				  {
				    long int gp_id=gp_ids[i];
				    long int gpp_id=gpp_ids[i];
				    scalar_type pp=p_part[i];
				    scalar_type S_pf_ppg=q[gp_id][t][f]*q[gpp_id][t][g]*pp;
				    scalar_type S_ppf_pg=q[gpp_id][t][f]*q[gp_id][t][g]*pp;
				    //S EVENT, events #1 and #2 in part c of Fig.A1 in http://arxiv.org/abs/1211.4606
				    //q[g_id][t][e]+=q[gp_id][t][f]*q[gpp_id][t][g] +q[gpp_id][t][f]*q[gp_id][t][g];
				    q_sum+= S_pf_ppg + S_ppf_pg;
				    //S.
				  }
			      q[g_id][t][e]=q_sum;

			    }

			  //branches that cross to next time slice
			  else
			    {
			      //trivial
			      ;//q[g_id][t][e]=q[g_id][t][e];
			    }
			}
		      //boundaries for branch e.
		    }
		}

	      if(1)
		{

		  //events within slice rank at time t on alpha virtual branch
		  scalar_type G_bar=Ge[-1][t];
		  //note that the coalescent approximation in http://arxiv.org/abs/1211.4606 is exp(-(Delta_bar*(n-N)/N+Lambda_bar)*Delta_t );

		  q[g_id][tpdt][alpha]=0;
		  scalar_type q_sum=0;
		  scalar_type q_sum_nl=0;
		  for (int branch_i=0;branch_i<n;branch_i++)
		    {
		      int e = time_slices[rank][branch_i];
		      scalar_type tau_e=vector_parameter["tau"][e];
		      scalar_type p_Ntau_e=tau_e*Delta_t;

		      //non-leaf directed partition
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {
			    long int gp_id=gp_ids[i];
			    long int gpp_id=gpp_ids[i];
			    scalar_type pp=p_part[i];
			    scalar_type T_ep_app=p_Ntau_e*q[gp_id][t][e]*q[gpp_id][t][alpha]*pp;
			    scalar_type T_ap_epp=p_Ntau_e*q[gp_id][t][alpha]*q[gpp_id][t][e]*pp;
			    //T EVENT, events #3 and #4 in part b of Fig.A1 in http://arxiv.org/abs/1211.4606
			    //q[g_id][tpdt][alpha]+=p_Ntau_e*(q[gp_id][t][e]*q[gpp_id][t][alpha]+q[gp_id][t][alpha]*q[gpp_id][t][e]);
			    q_sum_nl+=T_ep_app+T_ap_epp;
			    //T.

			  }
		    }
		  //non-leaf directed partition
		  if (not is_a_leaf)
		    for (int i=0;i<N_parts;i++)
		      {
			long int gp_id=gp_ids[i];
			long int gpp_id=gpp_ids[i];
			scalar_type pp=p_part[i];

			scalar_type Sb=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha])*pp;
			//S_bar EVENT, event #2 in part b of Fig.A1 in http://arxiv.org/abs/1211.4606
			//(note that Delta_bar corresponds to sigma, the Delta_bar,Lambda_bar distinction keeps track of speciaiton (birth) vs extiction (death),
			// but for the Moran process Delta_bar=Lambda_bar=sigma )
			//q[g_id][tpdt][alpha]+=p_Delta_bar*(2*q[gp_id][t][alpha]*q[gpp_id][t][alpha]);
			q_sum_nl+=Sb;
			//S_bar.

		      }

		  q[g_id][tpdt_nl][alpha]+=q_sum_nl;

		  for (int branch_i=0;branch_i<n;branch_i++)
		    {
		      int e = time_slices[rank][branch_i];
		      scalar_type tau_e=vector_parameter["tau"][e];
		      scalar_type p_Ntau_e=tau_e*Delta_t;
		      scalar_type TLb=p_Ntau_e*Ebar*q[g_id][t][e];
		      //TL_bar EVENT, event #5 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
		      //(note that since Ebar ~ 1, most transfers are expected to involve the TL evenet not the T event,
		      //this should not be confused with the TL event of the Tofigh/Doyon/ODTL models, which here corresponds
		      // to SL_bar + TL ..)
		      //q[g_id][tpdt][alpha]+=p_Ntau_e*Ebar*q[g_id][t][e];
		      q_sum+=TLb;
		      //TL_bar.

		    }
		  scalar_type empty=G_bar*q[g_id][t][alpha];
		  //0 EVENT, event #1 in part b of Fig.A1 in http://arxiv.org/abs/1211.4606
		  //q[g_id][tpdt][alpha]+=G_bar*q[g_id][t][alpha];
		  q_sum+=empty;
		  //0.

		  q[g_id][tpdt][alpha]+=q_sum;
		  //events within slice rank at time t on alpha virtual branch.
		}
	      if(1)
		{
 		  for (int branch_i=0;branch_i<n;branch_i++)
		    {
		      int e = time_slices[rank][branch_i];
		      scalar_type Get=Ge[e][t];
		      scalar_type Eet=Ee[e][t];
		      scalar_type delta_e=vector_parameter["delta"][e];
		      scalar_type p_delta_e=delta_e*Delta_t;

		      //events within slice rank at time t on branch e
		      q[g_id][tpdt][e]=0;
		      scalar_type q_sum=0;
		      scalar_type q_sum_nl=0;

		      //non-leaf directed partition
		      if (not is_a_leaf)
			for (int i=0;i<N_parts;i++)
			  {
			    long int gp_id=gp_ids[i];
			    long int gpp_id=gpp_ids[i];
			    scalar_type pp=p_part[i];
			    scalar_type qpe=q[gp_id][t][e];
			    scalar_type qppe=q[gpp_id][t][e];
			    scalar_type Sb_pa_ppe= p_Delta_bar*q[gp_id][t][alpha]*qppe*pp;
			    scalar_type Sb_pe_ppa= p_Delta_bar*qpe*q[gpp_id][t][alpha]*pp;

			    //S_bar EVENT, events #3 and #4 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
			    //(The majority of transfer events involve this event.)
			    //q[g_id][tpdt][e]+=p_Delta_bar*(q[gp_id][t][alpha]*q[gpp_id][t][e]+q[gp_id][t][e]*q[gpp_id][t][alpha]);
			    q_sum_nl+= Sb_pa_ppe + Sb_pe_ppa;
			    //S_bar.

			    scalar_type D=2*p_delta_e*qpe*qppe*pp;
			    //D EVENT, event #2 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
			    //q[g_id][tpdt][e]+=p_delta_e*q[gp_id][t][e]*q[gpp_id][t][e];
			    q_sum_nl+= D;
			    //D.

			  }

		      scalar_type SLb=p_Delta_bar*Eet*q[g_id][t][alpha];
		      //SL_bar EVENT, event #5 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
		      //(Transfer events where the donor copy is lost involve this event.)
		      //q[g_id][tpdt][e]+=p_Delta_bar*Eet*q[g_id][t][alpha];
		      q_sum_nl+=SLb;
		      //SL_bar.

		      q[g_id][tpdt_nl][e]+=q_sum_nl;
		      //if (q[g_id][tpdt_nl][e]>1) q[g_id][tpdt_nl][e]=1;//XX

		      scalar_type empty=Get*q[g_id][t][e];
		      //0 EVENT, event #1 in part a of Fig.A1 in http://arxiv.org/abs/1211.4606
		      //q[g_id][tpdt][e]=Get*q[g_id][t][e];
		      q_sum+=empty;
		      //0.

		      q[g_id][tpdt][e]+=q_sum;
                //if (q[g_id][tpdt][e]>1) q[g_id][tpdt][e]=1;
		      //events within slice rank at time t on branch e.
		    }
		}
	      //######################################################################################################################
	      //#########################################INNNER LOOP##################################################################
	      //######################################################################################################################
	    }
	}
      gp_ids.clear();
      gpp_ids.clear();
      p_part.clear();
    }

  scalar_type root_norm=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {
	      root_norm+=1;
	    }
	  root_norm+=1;
	}
    }

  scalar_type root_sum=0;
  for (int rank=0;rank<last_rank;rank++)
    {
      int n=time_slices[rank].size();
      for (int t_i=0;t_i<(int)time_slice_times[rank].size();t_i++)
	{
	  scalar_type t=time_slice_times[rank][t_i];
	  for (int branch_i=0;branch_i<n;branch_i++)
	    {
	      int e = time_slices[rank][branch_i];
	      //if (rank==last_rank-1 and t_i==(int)time_slice_times[rank].size()-1)//(1-Ee[e][time_slice_times[rank][t_i]])/
	      root_sum+=q[-1][t][e]/root_norm;
	    }
	  //if (rank==last_rank-1 and t_i==(int)time_slice_times[rank].size()-1)//(1-Ee[-1][time_slice_times[rank][t_i]]);
	  root_sum+=q[-1][t][alpha]/root_norm;
	}
    }

  //del-locs
  g_ids.clear();
  g_id_sizes.clear();

  return root_sum;
}

void exODT_model::calculate_EGb()
{

  for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ee.begin();it!=Ee.end();it++)//del_loc
    (*it).second.clear();
  Ee.clear();
  for (std::map<int,std::map <scalar_type,scalar_type> >::iterator it=Ge.begin();it!=Ge.end();it++)//del_loc
    (*it).second.clear();
  Ge.clear();


  map<int,scalar_type> Ee_y;//del-loc
  map<int,scalar_type> Ge_y;//del-loc
  map<int,scalar_type> E_k1,E_k2,E_k3,E_k4;//del-loc
  map<int,scalar_type> G_k1,G_k2,G_k3,G_k4;//del-loc

  map<int,double> tmp; //XX
  tmp[0]=1;
  tmp[1]=1;

  for (int rank=0;rank<last_rank;rank++)
    for (int tsi=0;tsi<(int)time_slice_times[rank].size();tsi++)
      {
	map<int,map <scalar_type,scalar_type> > y_E,y_G;//del-loc
	map<int,map <int,scalar_type> > iy_E,iy_G;//del-loc

	scalar_type t_b;
	if (tsi==(int)time_slice_times[rank].size()-1)
	  t_b = time_slice_begins[rank];
	else
	  t_b = time_slice_times[rank][tsi+1];
	scalar_type t_e;
	if (tsi==0)
	  {
	    if (rank>0 )
	      t_e = time_slice_begins[rank-1];
	    else
	      t_e = 0;
	  }
	else
	  {
	    t_e=time_slice_times[rank][tsi];
	  }
	scalar_type N=vector_parameter["N"][rank];

	scalar_type ni=time_slices[rank].size();
	scalar_type Delta_bar=vector_parameter["Delta_bar"][rank];//1
	scalar_type Lambda_bar=vector_parameter["Lambda_bar"][rank]*N/(N-ni);;
	scalar_type t=t_e;
	scalar_type tpdt=t_b;
	scalar_type h=(tpdt-t)/scalar_parameter["DD"];
	//scalar_type ti=t;
	scalar_type h_lambda_avg=h*scalar_parameter["lambda_avg"];
	scalar_type h_delta_avg=h*scalar_parameter["delta_avg"];
	scalar_type h_tau_avg=h*scalar_parameter["tau_avg"]*(N-ni)/(N-1)*N;
	scalar_type h_Delta_bar=h*Delta_bar;
	scalar_type h_Lambda_bar=h*Lambda_bar;


	for (int ii=0;ii<scalar_parameter["DD"];ii++)
	  {

	    //intial conditions
	    if (ii==0)
	      {
		if ( t==0)
		  Ee[-1][t]=1;
		//trivial else Ee[-1][t]=Ee[-1][t];

		//y_E[-1][t]=Ee[-1][t];
		iy_E[-1][ii]=Ee[-1][t];

		//Ee_y[-1]=y_E[-1][t];
		Ee_y[-1]=iy_E[-1][ii];

		Ge_y[-1]=1;

		//y_G[-1][t]=1;
		iy_G[-1][ii]=1;

	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		if (ii==0)
		  {
		    if ( t==0)
		      {
			Ee[e][t] = vector_parameter["fraction_missing"][e]; //0;
		      }
		    else if (t==t_end[e])
		      {
			int f=daughters[e][0];int g=daughters[e][1];
			Ee[e][t]=Ee[f][t]*Ee[g][t];
		      }
		    //trivial else{Ee[e][t]=Ee[e][t];}
		    //y_E[e][t]=Ee[e][t];
		    iy_E[e][ii]=Ee[e][t];

		    //Ee_y[e]=y_E[e][t];
		    Ee_y[e]=iy_E[e][ii];

		    Ge_y[e]=1;

		    //y_G[e][t]=1;
		    iy_G[e][ii]=1;

		  }
	      }
	    // RK4: 4th order Runge-Kutta for y'=f(y)
	    // k1 = f(y[n])
	    E_k1[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k1[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];
		E_k1[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k1[-1]-=h_tau_f*(1-Ee_y[f])* Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda;
		scalar_type h_delta=h*delta;

		// k1 = f(y[n])
		E_k1[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k1[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];

	      }
	    // k2 = f(y[n]+h/2 k1)

	    //Ee_y[-1]=y_E[-1][ti]+1/2.* E_k1[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1/2.* E_k1[-1];
	    //Ge_y[-1]=y_G[-1][ti]+1/2.* G_k1[-1];
	    Ge_y[-1]=iy_G[-1][ii]+1/2.* G_k1[-1];

	    E_k2[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k2[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];
		E_k2[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k2[-1]-=h_tau_f*(1-Ee_y[f])* Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda;
		scalar_type h_delta=h*delta;

		// k2 = f(y[n]+h/2 k1)
		//Ee_y[e] =y_E[e][ti]+1/2. * E_k1[e];
		Ee_y[e] =iy_E[e][ii]+1/2. * E_k1[e];
		//Ge_y[e] =y_G[e][ti]+1/2. * G_k1[e];
		Ge_y[e] =iy_G[e][ii]+1/2. * G_k1[e];

		E_k2[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k2[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];
	      }

	    // k3 = f(y[n]+h/2 k2)
	    //Ee_y[-1]=y_E[-1][ti]+1/2.* E_k2[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1/2.* E_k2[-1];
	    //Ge_y[-1]=y_G[-1][ti]+1/2.* G_k2[-1];
	    Ge_y[-1]=iy_G[-1][ii]+1/2.* G_k2[-1];

	    E_k3[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k3[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];
		E_k3[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k3[-1]-=h_tau_f*(1-Ee_y[f])* Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda;
		scalar_type h_delta=h*delta;

		// k3 = f(y[n]+h/2 k2)
		//Ee_y[e] =y_E[e][ti]+1/2. * E_k2[e];
		Ee_y[e] =iy_E[e][ii]+1/2. * E_k2[e];
		//Ge_y[e] =y_G[e][ti]+1/2. * G_k2[e];
		Ge_y[e] =iy_G[e][ii]+1/2. * G_k2[e];

		E_k3[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k3[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];
	      }

	    // k4 = f(y[n]+h k3)
	    //Ee_y[-1]=y_E[-1][ti]+1* E_k3[-1];
	    Ee_y[-1]=iy_E[-1][ii]+1* E_k3[-1];
	    //Ge_y[-1]=y_G[-1][ti]+1* G_k3[-1];
	    Ge_y[-1]=iy_G[-1][ii]+1* G_k3[-1];

	    E_k4[-1]=(h_Lambda_bar+h_lambda_avg)*(1-Ee_y[-1])-( h_Delta_bar+h_delta_avg+h_tau_avg)*(1- Ee_y[-1])* Ee_y[-1];
	    G_k4[-1]=-(( h_Delta_bar+h_delta_avg+h_tau_avg)*(1-2*Ee_y[-1]) + (h_Lambda_bar+h_lambda_avg) ) * Ge_y[-1];
	    for (int j=0;j<(int)time_slices[rank].size();j++)
	      {
		int f=time_slices[rank][j];
		scalar_type h_tau_f=h*vector_parameter["tau"][f];
		E_k4[-1]-=h_tau_f*(1-Ee_y[f])* Ee_y[-1];
		G_k4[-1]-=h_tau_f*(1-Ee_y[f]) *Ge_y[-1];
	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		scalar_type delta=vector_parameter["delta"][e];
		scalar_type lambda=vector_parameter["lambda"][e];
		scalar_type h_lambda=h*lambda;
		scalar_type h_delta=h*delta;

		// k4 = f(y[n]+h k3)
		//Ee_y[e] =y_E[e][ti]+1 * E_k3[e];
		Ee_y[e] =iy_E[e][ii]+1 * E_k3[e];

		//Ge_y[e] =y_G[e][ti]+1 * G_k3[e];
		Ge_y[e] =iy_G[e][ii]+1 * G_k3[e];

		E_k4[e]=h_lambda*(1-Ee_y[e])-( h_delta*(1- Ee_y[e]) + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1])) * Ee_y[e];
		G_k4[e]=-1*(h_lambda+h_delta*(1-2*Ee_y[e])  + (h_Delta_bar+h_tau_avg)*(1-Ee_y[-1]))* Ge_y[e];
	      }
	    // y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4)
	    //y_E[-1][ti+h]=Ee_y[-1] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    //iy_E[-1][ii+1]=Ee_y[-1] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);

	    ///*
	    if (ii==0)
	      iy_E[-1][ii+1]=Ee[-1][t] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    else
	      iy_E[-1][ii+1]=iy_E[-1][ii] + 1/6. * (E_k1[-1] + 2*E_k2[-1] + 2*E_k3[-1] + E_k4[-1]);
	    //*/
	    //y_G[-1][ti+h]=Ge_y[-1] + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);
	    //iy_G[-1][ii+1]=Ge_y[-1] + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);


	    if (ii==0)
	      iy_G[-1][ii+1]=1 + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);
	    else
	      iy_G[-1][ii+1]=iy_G[-1][ii] + 1/6. * (G_k1[-1] + 2*G_k2[-1] + 2*G_k3[-1] + G_k4[-1]);


	    if (ii==scalar_parameter["DD"]-1)
	      {
		//Ee[-1][tpdt]=y_E[-1][ti+h];
		Ee[-1][tpdt]=iy_E[-1][ii+1];

		//Ge[-1][t]=y_G[-1][ti+h];
		Ge[-1][t]=iy_G[-1][ii+1];

		//cout << -1 << " " << t << " " << Ee[-1][tpdt] << " " << Ge[-1][t]<<endl;

	      }

	    for (int i=0;i<(int)time_slices[rank].size();i++)
	      {
		int e=time_slices[rank][i];
		// y[n+1] = y[n] + h/6 (k1 + 2 k2 + 2 k3 + k4)
		//y_E[e][ti+h]=Ee_y[e] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		//iy_E[e][ii+1]=Ee_y[e] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);

		///*
		if (ii==0)
		  iy_E[e][ii+1]=Ee[e][t] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		else
		  iy_E[e][ii+1]=iy_E[e][ii] + 1/6. * (E_k1[e] + 2*E_k2[e] + 2*E_k3[e] + E_k4[e]);
		//*/

		//y_G[e][ti+h]=Ge_y[e] + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);
		//iy_G[e][ii+1]=Ge_y[e] + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);

		if (ii==0)
		  iy_G[e][ii+1]=1 + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);
		else
		  iy_G[e][ii+1]=iy_G[e][ii] + 1/6. * (G_k1[e] + 2*G_k2[e] + 2*G_k3[e] + G_k4[e]);

		if (ii==scalar_parameter["DD"]-1)
		  {
		    //Ee[e][tpdt]=y_E[e][ti+h];
		    Ee[e][tpdt]=iy_E[e][ii+1];
		    //Ge[e][t]=y_G[e][ti+h];
		    Ge[e][t]=iy_G[e][ii+1];

		    //if (e<2) cout << e << " " << t << " " << Ee[e][tpdt] << " " << Ge[e][t]<<" "<< tmp[e]<<endl;

		  }
	      }
	    //ti=ti+h;

	  }
      }
  //del-locs
  Ee_y.clear();
  Ge_y.clear();
  E_k1.clear();E_k2.clear();E_k3.clear();E_k4.clear();
  G_k1.clear();G_k2.clear();G_k3.clear();G_k4.clear();
  /*
  for (map<int,map <scalar_type,scalar_type> >::iterator it=y_E.begin();it!=y_E.end();it++)
    (*it).second.clear();
  y_E.clear();
  for (map<int,map <scalar_type,scalar_type> >::iterator it=y_G.begin();it!=y_G.end();it++)
    (*it).second.clear();
  y_G.clear();
  */
}
