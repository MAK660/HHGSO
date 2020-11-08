/*
 Team Formation Problem using Hybrid Henry Gas Optimizer
 Hybrid Algo: HGSO, Jaya Algorithm (JA), 
 Sooty Tern Optimization Algorithm (STOA), 
 Butterfly Optimization Algorithm (BOA) 
 and Owl Search Algorithm (OSA) 
 
 Author: Prof Kamal Z. Zamli
 Email: kamalz@ump.edu.my
 Creation Date: 14 Oct 2019 (12:02 am - Tamu Hotel, KL)
 The code is release as is. Please email me if you need help.
 Part of the codes are used in the following paper title:
 "Hybrid Henry Gas Solubility Optimization Algorithm with Dynamic Cluster-to-Algorithm Mapping"
 */
 

import java.io.*;
import java.util.*;
import java.text.DecimalFormat;


public class HHGSOv1
 {
 
     static String costs_file, skills_file, values_file;
   	 static ArrayList<String> skills_list = new ArrayList<String>();
	 static ArrayList<String> persons_list = new ArrayList<String>();
     static ArrayList<String> unprocess_skills_connections = new ArrayList<String>();
     static ArrayList<String> unprocess_costs_connections = new ArrayList<String>();
	 static ArrayList<String> memory_list = new ArrayList<String>();
	 static int total_persons;
	 static int total_skills;
	 static double[][] costs_connections;
     static int [][] skills_connections;
	 static int [] skills_to_find;
	 static int [] clone_skills_to_find;
	 static int max_fitness_evaluation=5000;
	 static int count_fitness_evaluation=0;
	 static int max_iteration = 200;
	 static int max_run=1;
	 
	 // Population array
	 static int[][] population;
	 static String[] population_short_seq_list; 
	 static int [] obj_value1_array; 	 
	 static double [] obj_value2_array;
	 static int population_size=30;
	 static String running_algorithm="";
	 
	
	 // Original Henry Gas variables are global
     static double I1=5E-2;
	 static double I2=100.0;
	 static double I3=1E-2;
	 static double T_theta = 298.15;
	 static int division = 5;
     static double henry_coefficient[] = new double [division];
	 static double partial_pressure[] = new double [division];
	 static double gas_constant[] = new double [division];
	 static int best_idx_value_in_cluster [] = new int [division]; 
	 static int best_idx_value_overall;
	 static double solubility; 
     static int N_size = (int)population_size/(division);
	 static double epsilon =0.05;
	 static double K_constant = 1.0;
	 static double alpha_constant = 1.0;
	 static double beta_constant = 1.0;
	 static double gamma_hg;
	 static double F_signed =1.0;		
	 
	
	 // Ensemble Hybrid Variables
	 static int division_=1;
	 static int N_size_ = (int) population_size/division_;
	 //static String defined_algo_list[] ={"Jaya","Cuckoo","SootyTern","Owl","Butterfly","HenryGas"};
	 static String defined_algo_list[] ={"Jaya","HenryGas","SootyTern","Owl","Butterfly"};
	 static int count_JA=0;
	 static int count_HGSO=0;
	 static int count_BOA=0;
     static int count_OSA=0;
     static int count_STOA=0;
	 
	 
	 static int best_idx_value_in_cluster_ [] = new int [division_]; 
	 static int best_idx_value_overall_;
	 static int poor_idx_value_in_cluster_ [] = new int [division_]; 
	 static int poor_idx_value_overall_;
	 
	 
   public static void main(String[] args) throws IOException
    {
	
		System.out.println ("TEAM FORMING OPTIMIZER ");
		process_cmd_line(args);
		
		load_skills_in_memory(skills_file);
		process_skills();
		initialize_skills_connections();	
		load_costs_in_memory(costs_file);
		process_costs();
		initialize_costs_connections();	
        load_skills_to_find(values_file);
		
		
		//-----------------------------------
		// Automatic Running Data Collection
		//-----------------------------------f
		double best_cost;
	    double array_of_best_cost[] = new double [max_run];
		String best_seq_str="";
		double std_dev;
		double mean_cost;
		double mean_elapsed_time;
		String array_of_best_team_members[] = new String [max_run];
		double worst_cost;
	    double array_of_worst_cost[]= new double [max_run];
		String worst_seq_str="";
		int total_member;
		int array_of_total_members[]=new int [max_run];
		long start_time;
        long end_time;;    
		double array_of_elapsed_time[]=new double[max_run];

        // measure effects of cluster size on mapping of meta-heuristics		
		
		double array_of_mean_count_JA[]=new double [max_run];
		double array_of_mean_count_HGSO[]=new double [max_run];
		double array_of_mean_count_BOA[]=new double [max_run];
		double array_of_mean_count_OSA[]=new double [max_run];
		double array_of_mean_count_STOA[]=new double [max_run];
		double mean_JA;
		double mean_HGSO;
	    double mean_BOA;
		double mean_OSA;
		double mean_STOA;
		  
		for (int run=0;run< max_run;run++) 
 		  {
             start_time = System.nanoTime();
			 
		 	 search_hybrid_henry_gas();
			 end_time = System.nanoTime();
			 double elapsed_time = (end_time-start_time)/1000000000.0;
			 
		     int x = get_index_best_cost();	
			 int y = get_index_worst_cost();
	         String best_seq_string= population_short_seq_list[x];
             total_member= obj_value1_array[x];
	         best_cost = obj_value2_array[x];
			 worst_cost = obj_value2_array[y];
			 System.out.println ("=======================================================");	  
	         System.out.println ("Run No = "+run);
			 System.out.println ("Execution time (sec) = "+elapsed_time);
		     System.out.println ("Algorithm = "+running_algorithm);
	         System.out.println ("Best Sequence ==> "+best_seq_string);
	         System.out.println ("No of persons/members   = "+total_member);
	         System.out.println ("Best cost value  = "+best_cost);
			 String S=remap_team_members(best_seq_string);
	         System.out.println (" [ Team Members ] "+S);
		     
			 array_of_best_cost[run]=best_cost;
			 array_of_worst_cost[run]=worst_cost;
			 array_of_best_team_members[run]=S;
		     array_of_total_members[run]=total_member;
			 array_of_elapsed_time [run]=elapsed_time;
		
		     int total_count = count_JA+count_HGSO+count_BOA+count_OSA+count_STOA;
			 
			 array_of_mean_count_JA[run]=(double)100*count_JA/total_count;
		     array_of_mean_count_HGSO[run]=(double)100*count_HGSO/total_count;
		     array_of_mean_count_BOA[run]=(double)100*count_BOA/total_count;
		     array_of_mean_count_OSA[run]=(double)100*count_OSA/total_count;
		     array_of_mean_count_STOA[run]=(double)100*count_STOA/total_count;
			 
             // reset count for each algorithm in hybrid henry gas 
			 count_JA=0;
		     count_HGSO=0;
		     count_BOA=0;
             count_OSA=0;
             count_STOA=0;			
			 
		  }

         for (int run=0;run< max_run;run++)
		  {
			 System.out.println ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");	  
             System.out.println ("Summary "+running_algorithm+ " Run = "+run);
		     System.out.println ("Execution time (sec) = "+array_of_elapsed_time[run]);
		     System.out.println ("Best Cost = "+array_of_best_cost[run]);
			 System.out.println ("Worst Cost = "+array_of_worst_cost[run]);
			 System.out.println ("No of Team members = "+array_of_total_members[run]);
			 System.out.println (" [ Team Members ] "+ array_of_best_team_members[run]);
			 System.out.println ("******Percentage Execution Algorithm Distribution******");
			 System.out.println ("Jaya Execution Distribution = "+array_of_mean_count_JA[run]);
		     System.out.println ("Henry Gas Execution Distribution = "+array_of_mean_count_HGSO[run]);
			 System.out.println ("Butterfly Execution Distribution = "+array_of_mean_count_BOA[run]);
			 System.out.println ("Owl Execution Distribution = "+array_of_mean_count_OSA[run]);
			 System.out.println ("Sooty Tern Execution Distribution = "+array_of_mean_count_STOA[run]);
		  }			 
 		  
		  
		  int index_best_of_the_best= 0;
		  double best_of_the_best=array_of_best_cost[0];
		  
		  for (int run=0;run<max_run;run++)
		    {
			   if (array_of_best_cost[run]<best_of_the_best)
			   {
				   index_best_of_the_best=run;
				   best_of_the_best=array_of_best_cost[run];
               }				   
			}   
			   
		  int index_worst_of_the_worst = 0;
		  double worst_of_the_worst=array_of_worst_cost[0];
		  
		  for (int run=0;run<max_run;run++)
		    {
			   if (array_of_worst_cost[run]>worst_of_the_worst)
			   {
				   index_worst_of_the_worst=run;
				   worst_of_the_worst=array_of_worst_cost[run];
               }				   
			}
		  
		  mean_cost=0.0;
		  for (int run=0;run<max_run;run++)
		    {
			    mean_cost=array_of_best_cost[run]+mean_cost;
			}
		  mean_cost=mean_cost/max_run;	
		  
		  std_dev=calculate_stdv(array_of_best_cost,mean_cost);
			
		  mean_elapsed_time=0.0;
		  for (int run=0;run<max_run;run++)
		    {
			    mean_elapsed_time=array_of_elapsed_time[run]+mean_elapsed_time;
			}
		  mean_elapsed_time=mean_elapsed_time/max_run;	
		 
		  mean_JA=0.0;
		  mean_HGSO=0.0;
		  mean_BOA=0.0;
		  mean_OSA=0.0;
		  mean_STOA=0.0;
		  for (int run=0;run<max_run;run++)
		    {
			    mean_JA=array_of_mean_count_JA[run]+mean_JA;
				mean_HGSO=array_of_mean_count_HGSO[run]+mean_HGSO;
				mean_BOA=array_of_mean_count_BOA[run]+mean_BOA;
				mean_OSA=array_of_mean_count_OSA[run]+mean_OSA;
				mean_STOA=array_of_mean_count_STOA[run]+mean_STOA;
			}
		  double total_mean= mean_JA+mean_HGSO+mean_BOA+mean_OSA+mean_STOA;
		  mean_JA=100*mean_JA/total_mean;	
		  mean_HGSO=100*mean_HGSO/total_mean;	
		  mean_BOA=100*mean_BOA/total_mean;	
		  mean_OSA=100*mean_OSA/total_mean;
          mean_STOA=100*mean_STOA/total_mean;

		  
          System.out.println ("_________________________________________________________________"); 		 
		  System.out.println ("=================================================================");	  
          System.out.println ("Final Overall Execution Summary => "+running_algorithm);
		  System.out.println ("Input File => "+values_file);
		  System.out.println ("Cluster (only for HHGSO) => "+division_);
		  System.out.println ("Population Size => "+population_size);
		  System.out.println ("Max Iteration => "+max_iteration);
		  System.out.println ("Max Fitness Evaluation => "+max_fitness_evaluation);
		  System.out.println ("_________________________________________________________________"); 		 
		  System.out.println ("=================================================================");	  
          System.out.println ("No of runs = "+max_run);
		  System.out.println ("Best Overall Cost = "+array_of_best_cost[index_best_of_the_best]);
		  System.out.println ("Worst Overall Cost = "+array_of_worst_cost[index_worst_of_the_worst]);
		  System.out.println ("Mean Overall Cost = "+mean_cost);
		  System.out.println ("Standard deviation = "+std_dev);
		  System.out.println ("Mean execution time (sec) = "+mean_elapsed_time);
		  System.out.println ("Best no of team members = "+array_of_total_members[index_best_of_the_best]);
		  System.out.println (" [ Best Team Members ] "+ array_of_best_team_members[index_best_of_the_best]);
		  System.out.println ("[******Percentage Mean Overall Execution Algorithm Distribution******]");
		  System.out.println ("Mean Henry Gas Execution Distribution = "+roundTwoDecimals(mean_HGSO));
		  System.out.println ("Mean Jaya Execution Distribution = "+roundTwoDecimals(mean_JA));
		  System.out.println ("Mean Sooty Tern Execution Distribution = "+roundTwoDecimals(mean_STOA));
		  System.out.println ("Mean Butterfly Execution Distribution = "+roundTwoDecimals(mean_BOA));
		  System.out.println ("Mean Owl Execution Distribution = "+roundTwoDecimals(mean_OSA));
		 
		
    }
	
	
	//////////////////////////////////////////////////////////////
    //   Initialize population 
    //////////////////////////////////////////////////////////////	
	public static void initialize_population()
	{
	    int arr[];
		int seq[];
		obj_value1_array = new int [population_size];
		obj_value2_array = new double [population_size];
		
	    population = new int [population_size][total_persons];
		population_short_seq_list = new String [population_size];
	    for (int row=0;row<population_size;row++)
		  {
		     arr=generate_random_sequence(total_persons);
			 seq = objective_function_(arr);
			
			 // long sequence - all
	         //System.out.println ("Long sequence => "+array_sequence_to_string(arr));
			 for (int col=0;col<total_persons;col++)
			   {
		           population[row][col]= arr[col];
               }
			   
		     //System.out.println ("Short sequence => "+array_sequence_to_string(seq));
			 // short seq - best only
	         population_short_seq_list[row]=array_sequence_to_string(seq);
			    
			 int obj_person=objective_value1_(seq);
		     double obj_costs=objective_value2_(seq); 
			 System.out.println ("--------------------------------------------------");
			 System.out.println ("Population no = "+row);
			 System.out.println ("Long sequence = "+array_sequence_to_string(arr));
			 System.out.println ("Short sequence = "+array_sequence_to_string(seq));
			 System.out.println ("No of person = "+obj_person);
             System.out.println ("Costs = "+obj_costs);
             obj_value1_array[row]=obj_person; 	 
	         obj_value2_array[row]=obj_costs;  
		  }
	
	}
	
	
    //////////////////////////////////////////////////////////////
    //   Initialize population
    //////////////////////////////////////////////////////////////	
	public static void initialize_population_henry_gas()
	{
	
	    int arr[];
		int seq[];
		Random r = new Random();
		obj_value1_array = new int [population_size];
		obj_value2_array = new double [population_size];
		
	    population = new int [population_size][total_persons];
		population_short_seq_list = new String [population_size];
		int idx=0;
	    for (int row=0;row<population_size;row++)
		  {
			 
			 // divide population in a cluster division
			 if (row==0)
              {
			    // Henry Gas 				
                best_idx_value_overall=row;
                best_idx_value_in_cluster[idx]=row;  
                henry_coefficient[idx] = I1*r.nextDouble();
	            partial_pressure[idx] = I2*r.nextDouble();
	            gas_constant[idx] = I3*r.nextDouble(); 
              }
			 else if (row%N_size==0)
              {  
		        
				// Henry Gas
				idx++;
				best_idx_value_in_cluster[idx]=row;  
                henry_coefficient[idx] = I1*r.nextDouble();
	            partial_pressure[idx] = I2*r.nextDouble();
	            gas_constant[idx] = I3*r.nextDouble(); 
              }
			  
		     arr=generate_random_sequence(total_persons);
			 seq = objective_function_(arr);
			
			 // long sequence - all
	         //System.out.println ("Long sequence => "+array_sequence_to_string(arr));
			 for (int col=0;col<total_persons;col++)
			   {
		           population[row][col]= arr[col];
               }
			   
		     //System.out.println ("Short sequence => "+array_sequence_to_string(seq));
			 //short seq - best only
	         population_short_seq_list[row]=array_sequence_to_string(seq);
			    
			 int obj_person=objective_value1_(seq);
		     double obj_costs=objective_value2_(seq); 
			 System.out.println ("--------------------------------------------------");
			 System.out.println ("Population no = "+row);
			 System.out.println ("Long sequence = "+array_sequence_to_string(arr));
			 System.out.println ("Short sequence = "+array_sequence_to_string(seq));
			 System.out.println ("Costs no of person = "+obj_person);
             System.out.println ("Costs on interaction = "+obj_costs);
             obj_value1_array[row]=obj_person; 	 
	         obj_value2_array[row]=obj_costs;  
 	
   	        // capture index of the best obj value in a cluster
			 if (obj_value2_array[row]<obj_value2_array[best_idx_value_in_cluster[idx]])
				  best_idx_value_in_cluster[idx]=row;
			 
			 // capture index of the best overall obj value
			 if (obj_value2_array[row]<obj_value2_array[best_idx_value_overall])
				  best_idx_value_overall=row;			 
			 
		  }
		  
	}

    //////////////////////////////////////////////////////////////
    //   Initialize population hybrid henry gas
    //////////////////////////////////////////////////////////////	
	public static void initialize_population_hybrid_henry_gas()
	{
	
	    int arr[];
		int seq[];
		Random r = new Random();
		obj_value1_array = new int [population_size];
		obj_value2_array = new double [population_size];
		
	    population = new int [population_size][total_persons];
		population_short_seq_list = new String [population_size];
		int idx=0;
	    for (int row=0;row<population_size;row++)
		  {
			 
			 // divide population in a cluster division
			 if (row==0)
              {
                best_idx_value_overall_=row;
                best_idx_value_in_cluster_[idx]=row;  
                poor_idx_value_overall_=row;
                poor_idx_value_in_cluster_[idx]=row;  
                 
              }
			 else if ((row%N_size_==0)&& (division_>1))
              {  
		        idx++;
				best_idx_value_in_cluster_[idx]=row;  
				poor_idx_value_in_cluster_[idx]=row;  
              }
			  
		     arr=generate_random_sequence(total_persons);
			 seq = objective_function_(arr);
			
			 // long sequence - all
	         //System.out.println ("Long sequence => "+array_sequence_to_string(arr));
			 for (int col=0;col<total_persons;col++)
			   {
		           population[row][col]= arr[col];
               }
			   
		     //System.out.println ("Short sequence => "+array_sequence_to_string(seq));
			 //short seq - best only
	         population_short_seq_list[row]=array_sequence_to_string(seq);
			    
			 int obj_person=objective_value1_(seq);
		     double obj_costs=objective_value2_(seq); 
			 System.out.println ("--------------------------------------------------");
			 System.out.println ("Population no = "+row);
			 System.out.println ("Long sequence = "+array_sequence_to_string(arr));
			 System.out.println ("Short sequence = "+array_sequence_to_string(seq));
			 System.out.println ("Costs no of person = "+obj_person);
             System.out.println ("Costs on interaction = "+obj_costs);
             obj_value1_array[row]=obj_person; 	 
	         obj_value2_array[row]=obj_costs;  
 	
   	        // capture index of the best min obj value in a cluster
			 if (obj_value2_array[row]<obj_value2_array[best_idx_value_in_cluster_[idx]])
				  best_idx_value_in_cluster_[idx]=row;
			 
   	        // capture index of the poor max obj value in a cluster
			 if (obj_value2_array[row]>obj_value2_array[poor_idx_value_in_cluster_[idx]])
				  poor_idx_value_in_cluster_[idx]=row;
			  		  
			 // capture index of the best min overall obj value
			 if (obj_value2_array[row]<obj_value2_array[best_idx_value_overall_])
				  best_idx_value_overall_=row;			 
			 
			 // capture index of the poor max overall obj value
			 if (obj_value2_array[row]>obj_value2_array[poor_idx_value_overall_])
				  poor_idx_value_overall_=row;
		  }
		  
	}
	
	
	//////////////////////////////////////////////////////////////
    //   Display current population 
    //////////////////////////////////////////////////////////////	
	public static void display_population()
	{
	
	  System.out.println ("-------------------DISPLAY POPULATION-----------------------------");
	  for (int i=0;i<population_size;i++)
	    {
	  
	      System.out.println ("Population no = "+i);
		  System.out.println ("Long sequence = "+array_sequence_to_string(get_population_long_sequence(i)));
		  System.out.println ("Short sequence = "+population_short_seq_list[i]);
		  System.out.println ("No of person = "+obj_value1_array[i]);
          System.out.println ("Costs = "+obj_value2_array[i]);
	   }
	
	
	}
	//////////////////////////////////////////////////////////////
    //   Get the ith population long sequence
    //////////////////////////////////////////////////////////////	
	public static int [] get_population_long_sequence (int idx)
	{
	   int arr[]=new int[total_persons];
	   
	   for (int col=0;col<total_persons;col++)
	      arr[col]=population[idx][col];
	
	   return (arr);
	}

   
	//////////////////////////////////////////////////////////////
    //   Get index for best minimum cost
    //////////////////////////////////////////////////////////////		
	public static int get_index_best_cost()
	{
	   
	   int idx=0;
	   double best_cost=Double.MAX_VALUE;
	   
	   for (int i=0;i<population_size;i++)
	      {
		     if (obj_value2_array[i]<best_cost)
			   {
			      best_cost=obj_value2_array[i];
			      idx=i;
			   }
		  }
	
	  return idx;
	}
	
	
	//////////////////////////////////////////////////////////////
    //   Get index for worst cost
    //////////////////////////////////////////////////////////////		
	public static int get_index_worst_cost()
	{
	   
	   int idx=0;
	   double worst_cost=Double.MIN_VALUE;
	   
	   for (int i=0;i<population_size;i++)
	      {
		     if (obj_value2_array[i]>worst_cost)
			   {
			      worst_cost=obj_value2_array[i];
			      idx=i;
			   }
		  
		  }
	
	  return idx;
	}
	
		
	//////////////////////////////////////////////////////////////
    //  Hybrid Henry Gas Solubility Optimization Algorithm
    //////////////////////////////////////////////////////////////	
	public static void search_hybrid_henry_gas()
	{
	   
	  int loop=0;
	  int seq[];
	  int current_seq[];
	  int pBest_long_seq[];
	  int updated_long_seq[]=new int[total_persons];
	  String updated_short_seq_string;
	  String current_short_seq_string;
 	  int updated_obj1;
      double updated_obj2;
      int current_obj1;
	  double current_obj2;
	  String current_running_algo=defined_algo_list[0];;
      double f_best_old,f_best_new;
	  double p_min=0.2, p_max=0.8;
	  
	  // for Hybrid Henry Gas only
	  double I1_=5E-2;
	  double I2_=100.0;
	  double I3_=1E-2;
	  double T_theta_ = 298.15;
	  double henry_coefficient_;
      double partial_pressure_ ;
	  double gas_constant_;
	  double solubility_; 
      double epsilon_ =0.05;
	  double K_constant_ = 1.0;
	  double alpha_constant_ = 1.0;
	  double beta_constant_ = 1.0;
	  double gamma_hg_;
	  double F_signed_ =1.0;	
	 
 	  Random r = new Random();
	  running_algorithm="Hybrid Henry Gas Solubility Optimization Algorithm";
	  
	  // at the start, randomize sequence of execution of algorithms
	  int mylist[]=generate_random_sequence (defined_algo_list.length);
   	  String clone_defined_algo_list[]=defined_algo_list.clone();
	  for (int c=0; c<mylist.length;c++)
	     {
		  defined_algo_list[c]=clone_defined_algo_list[mylist[c]];
		 }
		 
	  if (division_>defined_algo_list.length)
	    {
		   System.out.println ("Exiting program. Too many clusters of algorithm (i.e. insufficient algo definition)");
		   System.out.println ("No of defined algorithm = "+defined_algo_list.length);
		   System.out.println ("No of cluster = "+division);
		   System.exit(0);
	    }
	  initialize_population_hybrid_henry_gas();
	  
	  // for intialization purpose only
	  f_best_old=obj_value2_array[0];
	  f_best_new=obj_value2_array[0];
	  
	  while (loop<max_iteration)
 	   {
		  int idx=0;
		  
		  henry_coefficient_ = I1*r.nextDouble();
	      partial_pressure_ = I2*r.nextDouble();
	      gas_constant_ = I3*r.nextDouble(); 

	      for (int i=0;i<population_size;i++)   
			{
			   
			   // current population
		       current_seq=get_population_long_sequence(i);
		       current_short_seq_string=population_short_seq_list[i];
			   current_obj1=obj_value1_array[i];
	           current_obj2 = obj_value2_array[i];
			   System.out.println ("Cluster No => "+idx);
			   System.out.println ("Current seq => "+array_sequence_to_string(current_seq));
			   System.out.println ("Current short seq => "+current_short_seq_string);	
               System.out.println ("No of person [current] = "+current_obj1);
               System.out.println ("Costs [current] = "+current_obj2);
			   
			   if (i==0)
			    {
                    current_running_algo=defined_algo_list[idx];					
				    continue;
				}
			   else if((i%N_size_==0)&&(division_>1))
                 {
					idx++;
					current_running_algo=defined_algo_list[idx];					
                 } 
			
			
	           // Our defined algo:{"Jaya","Cuckoo","SootyTern","Owl","Butterfly"}with dynamic arrangements	
			   if(current_running_algo.equals("Jaya"))
			      {
					 count_JA++;
                  	 System.out.println ("Cluster = "+idx+" <==> Running algorithm = "+current_running_algo);
                     
			         int idx_best = best_idx_value_in_cluster_[idx];; 
	                 int best_long_seq[] = get_population_long_sequence(idx_best);
			   
			         int idx_worst = poor_idx_value_in_cluster_[idx];
			         int worst_long_seq[] = get_population_long_sequence(idx_worst);
			   
	                 for (int j=0;j<current_seq.length;j++)
				       {
					     updated_long_seq[j]=current_seq[j]+(int)(r.nextDouble()*(best_long_seq[j]-current_seq[j])) 
					                                 - (int)(r.nextDouble()*(worst_long_seq[j]-current_seq[j]));	

                         if (updated_long_seq[j]>total_persons)
					        updated_long_seq[j]=0;
					     if (updated_long_seq[j]<0)
					        updated_long_seq[j]=total_persons;
													 
				       } 					 
				  }  
				else if(current_running_algo.equals("Cuckoo"))
				  {
                  	 System.out.println ("Cluster =  "+idx+" <==> Running algorithm = "+current_running_algo);
					 
			         for (int j=0;j<current_seq.length;j++)
				       {
				         updated_long_seq[j]=current_seq[j]+(int)(1^Math.round(LevyFlight()));	
					
					     if (updated_long_seq[j]>total_persons)
					        updated_long_seq[j]=0;
					     if (updated_long_seq[j]<0)
					        updated_long_seq[j]=total_persons;
				       }
		  		  }
				else if(current_running_algo.equals("SootyTern"))
				  {
					 count_STOA++;
					 double C_f=2.0;
	                 double S_a=2.0; 
	                 double C_b=0.5*r.nextDouble();
	                 double C_st;
	                 double M_st;
	                 double D_st;
	  
                  	 System.out.println ("Cluster = "+idx+" <==> Running algorithm = "+current_running_algo);
					 double k=(i/population_size)*(2.0*Math.PI);
					 int idx_best = best_idx_value_in_cluster_[idx];; 
	                 int best_long_seq[] = get_population_long_sequence(idx_best);
			   
	                 for (int j=0;j<current_seq.length;j++)
			           {
					      M_st= C_b*(best_long_seq[j]-current_seq[j]);
					      C_st= S_a*current_seq[j];
					      D_st= C_st+ M_st;
					      double radius=Math.exp(k);
					      double x_prime=radius*Math.sin(k);
					      double y_prime=radius*Math.cos(k);
					      double z_prime=radius*k;
					      updated_long_seq[j]=(int)(D_st*(x_prime+y_prime+z_prime)*best_long_seq[j]);
					
				       }
		  		  }
				else if(current_running_algo.equals("Owl"))
				  {					
                     count_OSA++;			  
	                 double normalized_intensity_Ii;
	                 double intensity_change_Ic;
	                 double distance_info_R;
	                 double probability_pm=0.5;
	                 double alpha;
	                 double beta=1.9-(1.9*loop/max_iteration);
		 
                  	 System.out.println ("Cluster = "+idx+" <==> Running algorithm = "+current_running_algo);
					
                     int idx_best = best_idx_value_in_cluster_[idx]; 
	                 int best_long_seq[] = get_population_long_sequence(idx_best);
			   
			         int idx_worst = poor_idx_value_in_cluster_[idx];
					 double denominator;
					 if (Math.abs(obj_value2_array[idx_best]-obj_value2_array[idx_worst])<0.01)
					    denominator=0.01;
					 else
                        denominator=obj_value2_array[idx_best]-obj_value2_array[idx_worst];						 
					 
			         normalized_intensity_Ii=(double)(current_obj2-obj_value2_array[idx_worst])/denominator;
					 
			         System.out.println ("Normalized intensity value Ii = "+normalized_intensity_Ii);
			         distance_info_R = Math.sqrt (Math.pow(current_obj1-obj_value1_array[idx_best],2)
			                               +Math.pow(current_obj2-obj_value2_array[idx_best],2));
			         if (distance_info_R-0.01<0) // avoid zero
				         distance_info_R=0.01;  
			         System.out.println ("Distance info  = "+distance_info_R);		   
			         intensity_change_Ic= (normalized_intensity_Ii/Math.pow(distance_info_R,2)) + r.nextDouble();
			         if (intensity_change_Ic-0.01<0) // avoid zero
				         intensity_change_Ic=0.01;  
			         System.out.println ("Intensity change Ic ="+intensity_change_Ic);
	                 
					 if (r.nextDouble()< probability_pm)
			            {
				          for (int j=0;j<current_seq.length;j++)
					       {
					          alpha = r.nextDouble()/2; 
					          updated_long_seq[j]=current_seq[j]+(int)(beta*intensity_change_Ic*(Math.abs((alpha*best_long_seq[j])-current_seq[j])));       
				           }				    
			   	        }
			         else
                        {
                           for (int j=0;j<current_seq.length;j++)
					       {
			                  alpha = r.nextDouble()/2;   
					          updated_long_seq[j]=current_seq[j]-(int)(beta*intensity_change_Ic*(Math.abs((alpha*best_long_seq[j])-current_seq[j])));       
				           }
						}   
							   
		  		  }
				else if(current_running_algo.equals("Butterfly"))
				  {
					 count_BOA++;
	                 double fragrance[] =  new double [population_size];
	                 double c=0.01; // sensory modal
	                 double I=total_persons; //stimuli - problem dependent 
					 
					 System.out.println ("Cluster = "+idx+" <==> Running algorithm = "+current_running_algo);
                     for (int row=0;row<population_size;row++)
		                 {
			               double a=0.1+(0.2*row/population_size);   
                           //double a= 0.1*Math.exp(1.0098*(row/population_size));
			               fragrance[row]=c*Math.pow(I,a);
		                 }	
					 
                     int idx_best = best_idx_value_in_cluster_[idx]; 
	                 int best_long_seq[] = get_population_long_sequence(idx_best);
			   
			         if (r.nextDouble()<0.8)
			           { 
				         for (int j=0;j<current_seq.length;j++)
				            {
					          double r_val=r.nextDouble();				   
					          double r_square=r_val*r_val;
					          int saved=current_seq[j];
					           updated_long_seq[j]=current_seq[j]+(int)((r_square*best_long_seq[j]-current_seq[j])*fragrance[i]);
					   
					          if (current_seq[j]>total_persons)
					             current_seq[j]=0;
					          if (current_seq[j]<0)
					             current_seq[j]=total_persons;
					  
					        }
				       }
					 else
					   { 			
                           			  
				         for (int j=0;j<current_seq.length;j++)
				           {
					          double r_val=r.nextDouble();	
				   	          double r_square=r_val*r_val;
					   
					          // find one peer butterfly
							  int k=0;
			                  while (true)
			                   {  
						          if (idx==0)
									  k=0;
								  else
									  k=idx*N_size_;
				                  int no=r.nextInt(N_size_);
								  k= k+no;
				                  if (k!=i)
					              break;
				               }
					   
                              int kth_seq[] = get_population_long_sequence(k); 
					          int saved=current_seq[j];
					          updated_long_seq[j]=current_seq[j]+(int)((r_square*kth_seq[j]-current_seq[j])*fragrance[i]);
					   
					          if (current_seq[j]>total_persons)
					             current_seq[j]=0;
					          if (current_seq[j]<0)
					             current_seq[j]=total_persons;  
                           }
                     							   
				       }	  
			      }
				else if(current_running_algo.equals("HenryGas"))
				  {
					  count_HGSO++;	   
				      double temperature_ = Math.exp(loop/max_iteration);
				      double val_ = -gas_constant_*(1/temperature_ - 1/T_theta_);
                      henry_coefficient_ = henry_coefficient_*Math.exp(val_);
				      solubility = K_constant_*henry_coefficient_* partial_pressure_; 
					  
					   
					  System.out.println ("Cluster = "+idx+" <==> Running algorithm = "+current_running_algo);
                      int idx1=best_idx_value_in_cluster_[idx];
			          int idx2=best_idx_value_overall_;
			   		
					
	                  for (int j=0;j<current_seq.length;j++)
				        {
				           if (r.nextDouble()>0.5)
					          F_signed=-1.0;
				           else
					          F_signed=1.0;
				  
				           gamma_hg = beta_constant*Math.exp(-(obj_value2_array[best_idx_value_overall]+epsilon)/(current_obj2+epsilon));
					       updated_long_seq[j]=current_seq[j]
					                    +(int)(F_signed*r.nextDouble()*gamma_hg*(population[idx1][j]-current_seq[j]))
					                    +(int)(F_signed*r.nextDouble()*alpha_constant*(solubility*population[idx2][j]-current_seq[j]));										
				        }
					  

				  }  
	            // must ensure legally correct sequence
                updated_long_seq = ensure_legal_sequence(updated_long_seq);
			    seq = objective_function_(updated_long_seq);
			  	updated_short_seq_string=array_sequence_to_string(seq);
		        updated_obj1=objective_value1_(seq);
		        updated_obj2=objective_value2_(seq);	
			    System.out.println ("Updated current seq => "+array_sequence_to_string(updated_long_seq));
			    System.out.println ("Updated short seq => "+updated_short_seq_string);	
                System.out.println ("No of person [updated] = "+updated_obj1);
                System.out.println ("Costs [updated] = "+updated_obj2);	
								
				
                // if the cost is improved (i.e. much less)  
                // replace the current population 
			    if (updated_obj2<current_obj2)
				  {
				     System.out.println ("Best Seq updated (i.e. updated<current)");
                     					                 
				     // update current population long sequence
					 for (int col=0;col<total_persons;col++)
				       population[i][col]=updated_long_seq[col];
						  
					 // update current population short sequence
					 population_short_seq_list[i] = updated_short_seq_string;
					
					// update objective values
					 obj_value1_array[i]=updated_obj1;
	                 obj_value2_array[i]=updated_obj2;
				  }
				 
   	            // capture index of the best min obj value in a cluster
			    if (obj_value2_array[i]<obj_value2_array[best_idx_value_in_cluster_[idx]])
				     best_idx_value_in_cluster_[idx]=i;
			 
   	            // capture index of the poor max obj value in a cluster
			    if (obj_value2_array[i]>obj_value2_array[poor_idx_value_in_cluster_[idx]])
				     poor_idx_value_in_cluster_[idx]=i;
			  		  
			    // capture index of the best min overall obj value
			    if (obj_value2_array[i]<obj_value2_array[best_idx_value_overall_])
				     best_idx_value_overall_=i;			 
			 
			    // capture index of the poor max overall obj value
			    if (obj_value2_array[i]>obj_value2_array[poor_idx_value_overall_])
				     poor_idx_value_overall_=i;
				
				
			}
            			
			// update Nx poor solution
	/*		int Nx = (int)(0.2*population_size);
		    for (int xx=0;xx<Nx;xx++)
 			   {
                 int idx_worst = get_index_worst_cost();
			     int worst_long_seq[] = get_population_long_sequence(idx_worst);
                 current_obj1=obj_value1_array[idx_worst];
				 current_obj2=obj_value2_array[idx_worst];
				 
				 // generate random sequence
				 updated_long_seq =generate_random_sequence(total_persons);
			     seq = objective_function_(updated_long_seq);
				 updated_short_seq_string=array_sequence_to_string(seq);
                 updated_obj2=objective_value2_(seq);	
				 updated_obj1=objective_value1_(seq);
			     if (updated_obj2<current_obj2)
				  {
				     System.out.println ("Poor Seq updated (i.e. updated<current)");
                     					                 
				     // update current population long sequence
					 for (int col=0;col<total_persons;col++)
				       population[idx_worst][col]=updated_long_seq[col];
						  
					 // update current population short sequence
					 population_short_seq_list[idx_worst] = updated_short_seq_string;
					
					// update objective values
					 obj_value1_array[idx_worst]=updated_obj1;
	                 obj_value2_array[idx_worst]=updated_obj2;
				  }
				 
		       }  */
			
			best_idx_value_overall_=get_index_best_cost();
		    poor_idx_value_overall_=get_index_worst_cost();
			f_best_new=obj_value2_array[best_idx_value_overall_];
			  
			
			// adaptive dynamic mapping of algo to clusters 
            if (loop==0)
  		   	   { 
				 f_best_old=obj_value2_array[best_idx_value_overall_];
				 f_best_new=obj_value2_array[best_idx_value_overall_];
			   }	
			else
			   {
                 if (f_best_new>=f_best_old) // non-improving co-operation 
				   {
					  System.out.println ("Non-Improving Moves - may change/maintain cluster-algo mapping with adaptive probability");
              
			          // Simple Adaptive Probability
            		  double p_threshold = p_max + (loop*(p_min-p_max)/max_iteration);
					  //double p_threshold = 1.0 - p_min + (loop*(p_max-p_min)/max_iteration);
					
					  System.out.println ("Probability threshold = "+p_threshold);
                      if (r.nextDouble()<p_threshold)   // probability whether to go to change cluster-algo mapping
                        {
                          System.out.println ("Must change current cluster algo mapping");
					      mylist=generate_random_sequence (defined_algo_list.length);
   	                      clone_defined_algo_list=defined_algo_list.clone();
						  for (int c=0; c<mylist.length;c++)
						    {
							  defined_algo_list[c]=clone_defined_algo_list[mylist[c]];
						    }
                         
                        }

                   } 	
				   
               }
			   
		    f_best_old=f_best_new;
			loop++;
		    if (count_fitness_evaluation>=max_fitness_evaluation)
		      break;
	   }
       
	    int idx = get_index_best_cost();	
	    String best_seq_string= population_short_seq_list[idx];
        int no_of_person = obj_value1_array[idx];
	    double best_cost = obj_value2_array[idx];
		
	    System.out.println ("--------------------------------------------");
	    System.out.println (running_algorithm);
	    System.out.println ("Final Best Sequence ==> "+best_seq_string);
	    System.out.println ("No of persons         = "+no_of_person);
	    System.out.println ("Obj Best Value [cost] = "+best_cost);
	    System.out.println (" [ Team Members ] "+remap_team_members(best_seq_string));		   
		
	}
	
	
	//////////////////////////////////////////////////////////////
    // Ensure legal sequence within the valid range/no repetition
    //////////////////////////////////////////////////////////////
    public static int [] ensure_legal_sequence(int arr[])
	{
	  int size=arr.length;
	  ArrayList<Integer> temp_list = new ArrayList<Integer>();
    
	  for (int i=0;i<size;i++)
	    {
		  
		  if (arr[i]<0)
		    arr[i]=Math.abs(arr[i]);
		  
		  if (arr[i]>size-1)
		    arr[i]=arr[i]%size;
			
		  if (!temp_list.contains(arr[i]))
		    temp_list.add(arr[i]);
		}
		
	  int ref[] = generate_random_sequence(size);
	  for(int i=0;i<size;i++)
	   {	      	
          if (!temp_list.contains(ref[i]))
		    temp_list.add(ref[i]);
	      
		  if (temp_list.size()==size)
		    break;	  
	   }
 	
	  int return_arr[]=new int[temp_list.size()];
	  for(int i=0;i<size;i++)
	   {
	     return_arr[i]=temp_list.get(i);
	   }
	   
	  return return_arr;
	  
	}
	
    //////////////////////////////////////////////////////////////
    // Synchronous add of each element in the sequence with val
    //////////////////////////////////////////////////////////////
    public static int [] add_sequence(int arr[], int val)
    {
        
		 int length = arr.length;
		 
         if (val>=length)
             val=val%length;
	     
		 // must force change at least 1 
		 if (val==0)
             val=1;
			 
         for (int i=0;i<length;i++)
            {
              arr[i]=arr[i]+val;
              if (arr[i]>=length)
                 arr[i]=arr[i]-length;                
            }

         return arr;

    }
	
	//////////////////////////////////////////////////////////////
    // Synchronous substract of each element in the sequence with val
    //////////////////////////////////////////////////////////////
    public static int [] substract_sequence(int arr[], int val)
    {
        
		 int length = arr.length;
		 
         if (val>=length)
             val=val%length;

		 // must force change at least 1 
		 if (val==0)
             val=1;
			 
         for (int i=0;i<length;i++)
            {
              arr[i]=arr[i]-val;
              if (arr[i]<0)
                 arr[i]=arr[i]+length;                
            }

         return arr;

    }
	
	//////////////////////////////////////////////////////////////
    //   Process command line 
    /////////////////////////////////////////////////////////////
	public static void process_cmd_line(String[] args) throws IOException
      {

        // Process argument lists  

        if (args.length == 0)
         {
            System.out.println("Missing problem definition files");
            System.exit(0);
         }
        else
         {
            for (int i=0;i< args.length;i++)
             {
              if (args[i].equals("-c")) 
                {
                    if (i + 1 < args.length) 
                     {
                        i++;
                        costs_file = args[i];
                     }
                  
                }
		       else if (args[i].equals("-s"))
                {
				    if (i + 1 < args.length) 
                     {
                        i++;
                        skills_file = args[i];
                     }
		        } 
				else if (args[i].equals("-v"))
                {
				    if (i + 1 < args.length) 
                     {
                        i++;
                        values_file = args[i];
                     }
		        }
				else if (args[i].equals("-r"))
                {
				    if (i + 1 < args.length) 
                     {
                        i++;
                        max_run = Integer.parseInt(args[i]);
                     }
		        }
                   
			  }	
          
         }
	 
	  }
	  
	////////////////////////////////////////////////////////////
    //   Display skills_connections metric
    ////////////////////////////////////////////////////////////
    public static void display_skills_connections()
	 {
	 
	    for (int i=0;i<total_persons;i++)
		   for (int j=0;j<total_skills;j++)
		      {
			      if (j<total_skills-1)
		          {
			        System.out.print(skills_connections[i][j]+":");
           	      } 
                 else
		          {
                    System.out.println(skills_connections[i][j]);
                  } 
			     
			  }
	 }
    
	  
	////////////////////////////////////////////////////////////
    //   Display costs_connections metric
    ////////////////////////////////////////////////////////////
    public static void display_costs_connections()
	 {
	 
	    for (int i=0;i<total_persons;i++)
		   for (int j=0;j<total_persons;j++)
		      {
			      if (j<total_persons-1)
		          {
			        System.out.print(costs_connections[i][j]+":");
           	      } 
                 else
		          {
                    System.out.println(costs_connections[i][j]);
                  } 
			     
			  }
	 }

	
	////////////////////////////////////////////////////////////
    //   Load skills to find values_file as array of int 
    ////////////////////////////////////////////////////////////
    public static void load_skills_to_find (String values_file) 
	                                            throws IOException
	{
	     RandomAccessFile f = new RandomAccessFile(values_file,"rw");
         long length = f.length();
         long position = 0;
         String content;
		 f.seek(0);
         ArrayList<Integer> match_idx = new ArrayList<Integer>();

	     while (position < length)
          {
		    content = f.readLine();
			String results[] = content.split("=");
			String val[] = results[1].split(",");
            for (int i=0;i<val.length;i++)
               {			
			     int index = skills_list.indexOf(val[i].trim());
				 System.out.println ("Skills to search for ==> "+val[i]);
				 if (index==-1)
				   {
                     System.out.println ("One of the skills requested does not exist...");
           			 System.exit(0);
				   } 				  
				 else
				   {
                      match_idx.add(index);
				   }
  
			   }	 
			position = f.getFilePointer();
		  
          }
        f.close();		 
		
		
		Collections.sort(match_idx); 
		skills_to_find= new int[total_skills];
		Arrays.fill(skills_to_find, 0);
		
		for (int i=0;i<total_skills;i++)
		   if (match_idx.contains(i))
			   skills_to_find[i]=1;   
		  
	}	 
	  
    ////////////////////////////////////////////////////////////
    //   Load skills in memory  
    ////////////////////////////////////////////////////////////
    public static void load_skills_in_memory (String skills_file) 
	                                            throws IOException
	{
	     RandomAccessFile f = new RandomAccessFile(skills_file,"rw");
         long length = f.length();
         long position = 0;
         String content;
		 int count=0;
	     System.out.println ("Loading skills file in memory => "+skills_file);
        // rewind file to position 0
         f.seek(0);
		 memory_list.clear();
         while (position < length)
          {
		    if (count%50==0)
		      System.out.println ("Count = ["+count+"] Loading skills to memory ...please wait");
            count++;
			content = f.readLine();
			memory_list.add(content);
	        position = f.getFilePointer();
		  }	
		 f.close(); 
	}
	
	////////////////////////////////////////////////////////////
    //   Load costs in memory  
    ////////////////////////////////////////////////////////////
    public static void load_costs_in_memory (String costs_file) 
	                                            throws IOException
	{
	     RandomAccessFile f = new RandomAccessFile(costs_file,"rw");
         long length = f.length();
         long position = 0;
         String content;
		 int count=0;
	     System.out.println ("Loading costs file in memory => "+skills_file);
        // rewind file to position 0
         f.seek(0);
		 memory_list.clear();
         while (position < length)
          {
		    if (count%50==0)
		      System.out.println ("Count = ["+count+"] Loading costs to memory ...please wait");
            count++;
			content = f.readLine();
			memory_list.add(content);
	        position = f.getFilePointer();
		  }	
		 f.close(); 
	}
	
    ////////////////////////////////////////////////////////////
    //   Process skills from memory 
    ////////////////////////////////////////////////////////////
    public static void process_skills () 
	{
	     String content;
	     System.out.println ("Processing skills ... ");
      
         for(int idx=0;idx<memory_list.size();idx++)
          {
            content = memory_list.get(idx);
			boolean process_already=false;
			// parse all the data for processing
            String results[] = content.split("=");
			for (int i=0;i<results.length;i++)
			   {
			      String string_val = results[i];
                  if (string_val.trim().equals("Total Persons".trim()))
				   {  
				       total_persons=Integer.parseInt(results[1]);
					   process_already=true;
				   }
			      else if (string_val.trim().equals("Total Skills".trim()))
				   {  
				       total_skills=Integer.parseInt(results[1]);
					   process_already=true;
				   }
			      else if (string_val.trim().equals("Persons Mapping".trim()))
                   {
                      String results2[]=results[1].split(",");
					  for (int j=0;j<results2.length;j++)
					    {
						  persons_list.add(results2[j]);
						  System.out.println ("Adding person => "+results2[j]);
						}
				      process_already=true;		
                   }      
				  else if (string_val.trim().equals("Skills Mapping".trim()))
                   {
                       String results2[]=results[1].split(",");
					   for (int j=0;j<results2.length;j++)
					    {
						  skills_list.add(results2[j]); 
						  System.out.println ("Adding skill => "+results2[j]);
						}
					   process_already=true;
                   } 
				  else // store unprocess skills connections 
				   {
				   
				       if (!unprocess_skills_connections.contains(content)&& process_already==false)
						    unprocess_skills_connections.add(content);
			       }
			   }
		  }	
	}
	
	////////////////////////////////////////////////////////////
    //   Process costs from memory 
    ////////////////////////////////////////////////////////////
    public static void process_costs () 
	{
	     String content;
	     System.out.println ("Processing costs ... ");
      
         for(int idx=0;idx<memory_list.size();idx++)
          {
            content = memory_list.get(idx);
			boolean process_already=false;
			// parse all the data for processing
            String results[] = content.split("=");
			for (int i=0;i<results.length;i++)
			   {
			      String string_val = results[i];
                  if (string_val.trim().equals("Persons Mapping".trim()))
                   {
				      process_already=true;
                      continue;	
                   }       
				  else // store unprocess costs connections 
				   {
				       if (!unprocess_costs_connections.contains(content)&& process_already==false)
 			            unprocess_costs_connections.add(content);
			       }
			   }
		  }	
	}
	
	//////////////////////////////////////////////////////////////
    //   Initialize skills connections 
    //////////////////////////////////////////////////////////////
     public static void initialize_skills_connections()
	  {
	  
		skills_connections = new int [total_persons][total_skills];
	     // initialize initial connections
		System.out.println ("Initializing skills connection metrics..");
        for (int row=0;row<total_persons;row++)
         for (int col=0;col<total_skills;col++)
           skills_connections[row][col]=0; 
   
    //    System.out.println ("Unprocess connections from skills file => "+unprocess_skills_connections.size());
		  
		for (int i=0;i<unprocess_skills_connections.size();i++)
           {
              String s = unprocess_skills_connections.get(i).trim();
			  String result[] = s.split("=");
			  String val[] = result[1].split(":");
	    	  for (int j=0;j<val.length;j++)
		        skills_connections[i][j]= (int)Integer.parseInt(val[j]);
           }    
		   
      }
	  
		
    ///////////////////////////////////////////////////////////
    //   Initialize cost metric 
    ///////////////////////////////////////////////////////////
     public static void initialize_costs_connections()
	  {
	  
	    costs_connections = new double [total_persons][total_persons];
	     // initialize initial connections
		System.out.println ("Initializing costs connection metrics..");
        for (int row=0;row<total_persons;row++)
         for (int col=0;col<total_persons;col++)
           costs_connections[row][col]=0.0; 
	    
		//System.out.println ("Unprocess connections from costs file => "+unprocess_costs_connections.size());
		
		   
        for (int i=0;i<unprocess_costs_connections.size();i++)
           {
              String s = unprocess_costs_connections.get(i).trim();
			  String result[] = s.split("=");
			  String val[] = result[1].split(":");
	    	  for (int j=0;j<val.length;j++)
		        costs_connections[i][j]= (double)Double.parseDouble(val[j]);
           } 
	  }


    //////////////////////////////////////////////////////////////
    // Calculate  objective value in terms the number of persons
    // only meaningful after objective function call
    /////////////////////////////////////////////////////////////
    public static int objective_value1_ (int seq[])
    {

      return (seq.length);

    }
	
    //////////////////////////////////////////////////////////////
    // Calculate  objective value in terms the connection costs
    // only meaningful after objective function call
    /////////////////////////////////////////////////////////////
    public static double objective_value2_ (int seq[])
    {
       double cost=0.0;
	 
	   for (int i=0;i<seq.length;i++)
	     {
		   int row = seq[i];
		   for (int j=i+1;j<seq.length;j++)
		     {
			   int col =seq[j];
			   cost = cost+costs_connections[row][col];
			 }
		 }
	   
	  return(cost); 
    }

	
    //////////////////////////////////////////////////////////////
    // Objective function - based on reduce sequence
    // concatenate long sequences of skills
    /////////////////////////////////////////////////////////////
    public static int [] objective_function_ (int arr [])
    {
	   boolean consist_skills_value=false;
	   int update_answer[] = new int[total_skills];
	   Arrays.fill(update_answer, 0);
	   
	   ArrayList<Integer> trim_seq = new ArrayList<Integer>();
	   String value_skills = array_sequence_to_string (skills_to_find);
		
	   clone_skills_to_find=skills_to_find.clone();
	
       for (int i=0;i<arr.length;i++)
        {
		  consist_skills_value=false;
          if (i==0)
           {
		      for (int j=0;j<total_skills;j++)
			     update_answer[j]=skills_connections[arr[i]][j];
				 
			  //System.out.println ("Current seq value = "+arr[i]);
			  //String s =array_sequence_to_string (update_answer);
              //System.out.println ("Update s       = "+s); 
			  //System.out.println ("Skills to find = "+value_skills);
              consist_skills_value=contain_skills_to_find(update_answer,clone_skills_to_find);
              if (consist_skills_value)   
			     trim_seq.add(arr[i]);
           }
          else
           {
              int tmp_answer[] = new int[total_skills];
	          Arrays.fill(tmp_answer, 0);
		      for (int j=0;j<total_skills;j++)
		         tmp_answer[j]=skills_connections[arr[i]][j];
			  
              //System.out.println ("Current seq value = "+arr[i]);
			  //String s =array_sequence_to_string (tmp_answer);
              //System.out.println ("Update s       = "+s); 
			  //System.out.println ("Skills to find = "+value_skills); 
			 	 
			  
              consist_skills_value=contain_skills_to_find(tmp_answer,clone_skills_to_find);
              if (consist_skills_value)
			    {
			      trim_seq.add(arr[i]);	
                  update_answer= merge_elements_(update_answer,tmp_answer);	   
				  if (complete_merge_sequence_(update_answer,skills_to_find))
		              break;	 
				} 
			    //s =array_sequence_to_string (update_answer);
        	    //System.out.println ("Skills so far  = "+s);  
             
           }
        } 
		
	  count_fitness_evaluation ++;	
	  int trim_result[] = new int [trim_seq.size()];
      for (int i=0;i<trim_seq.size();i++)
         trim_result[i]=trim_seq.get(i);
	
      return trim_result;
    }

    //////////////////////////////////////////////////////////
    //  Check if the given sequence has at least 1 skills 
    //  that we are looking for
    //////////////////////////////////////////////////////////
    public static boolean contain_skills_to_find (int arr[], int clone_skills_to_find[])
    {
      boolean outcome = false;
      for (int x=0; x<clone_skills_to_find.length; x++)
      {
       if (clone_skills_to_find[x]==1 && arr[x]==1)
		 {
		   outcome=true;
		   clone_skills_to_find[x]=0;  // to signify already covered skill
		   //break;
		 }
      }

      return outcome;
    }

      
	//////////////////////////////////////////////////////////
    //  Check for complete merge against skills to find
    //////////////////////////////////////////////////////////
    public static boolean complete_merge_sequence_ (int arr[], int skills_to_find[])
    {
      boolean outcome = true;
      for (int x=0; x<skills_to_find.length; x++)
      {
       if (skills_to_find[x]==1 && arr[x]==0)
		 {
		   outcome=false;
		   break;
		 }
      }

      return outcome;
    }
 

 
    ////////////////////////////////////////////////////////////
    //   Remap team members 
    ////////////////////////////////////////////////////////////
    public static String remap_team_members(String sequence)
	{
	   String[] list = sequence.split("-"); //split sequence
	   String team="\n";
       for (int i=0;i<list.length;i++)
         {
           int idx = Integer.parseInt(list[i]);
           String name = "===>"+persons_list.get(idx);
           if (i<list.length-1)
		     team=team+name+"\n";
		   else
             team=team+name;
		    
         } 
	  return team;
	}

    ////////////////////////////////////////////////////////////
    //   Generate random sequence 
    ////////////////////////////////////////////////////////////
    public static int [] generate_random_sequence (int dim)
    {
        int[] ar = new int[dim];
        int d, tmp;
        Random generator = new Random();
        int[] dummy_solution_array = new int[dim];
  
        for (int counter=0; counter<dim;counter++)
         {  
           dummy_solution_array [counter] = counter;
         }
       
        // copy array from seq_arr
        ar = dummy_solution_array.clone();

        // swap with new ar with random index
        // the first index and last index must not change
        for (int i=0;i<dim-1;i++)
        {
          d=i+(generator.nextInt()&(dim-1-i));
          tmp=ar[i];
          ar[i]=ar[d];
          ar[d]=tmp;
        }
       return ar;
    }
    
	 
	//////////////////////////////////////////////////////////////
    // Swap between two specified array points
    //////////////////////////////////////////////////////////////
    public static int [] swap_two_array_points (int point1, int point2,int arr[])
     {

         int t_arr[] = new int[arr.length];
         
         t_arr=arr.clone();
        
         int temp = t_arr[point1];
         t_arr[point1] = t_arr[point2];
         t_arr[point2] = temp;

         return t_arr;
     }


   	 
	//////////////////////////////////////////////////////////////
    //   k_opt array processing based on swap 2 points
    //////////////////////////////////////////////////////////////
    public static void k_opt_array_processing (int arr[])
     {

         int t_arr[] = new int[arr.length];      
        // t_arr=arr.clone();
		 
         for (int i=0;i<arr.length;i++)
		   {
		      for (int j=i+1;j<arr.length;j++)
			    {
				   t_arr=swap_two_array_points(i,j,arr);
				   //String s =array_sequence_to_string (t_arr);
				   //System.out.println ("K-Opt => "+s);
				}
		   }
     }


    /////////////////////////////////////////////////////////////
    //   Convert array sequence to string
    /////////////////////////////////////////////////////////////
   public static String array_sequence_to_string (int array [])
    {
       String sequence="";
       int length = array.length;
    
        for (int i=0;i<array.length;i++)
        {
         if (i<array.length-1)
            sequence=sequence+Integer.toString(array[i])+"-";
          else
            sequence=sequence+Integer.toString(array[i]);
        }

        return (sequence);
    }

   //////////////////////////////////////////////////////////////
   //   Convert string sequence to array
   /////////////////////////////////////////////////////////////
   public static int [] string_sequence_to_array (String sequence)
   {

     String result[] = sequence.split("-");
     int arr[] =new int [result.length];
     for (int i=0;i<result.length;i++)
	   {
		 arr[i]=Integer.parseInt(result[i]);
	   }	 
     return arr;
   }
  
   ////////////////////////////////////////////////////////////
   //    Merge elements
   ///////////////////////////////////////////////////////////
   public static int [] merge_elements_ (int arr1[],int arr2[])
   {
    
      int merge_array[] = new int[total_skills];
      Arrays.fill(merge_array, 0);
	
	 
     for (int x=0; x<arr1.length; x++)
     {

        if (arr1[x]==arr2[x])
          merge_array[x]=arr1[x];
        else if (arr1[x]==0 && arr2[x]==1)
          merge_array[x]=1;
        else if (arr1[x]==1 && arr2[x]==0)
          merge_array[x]=1;
     }
	 
     return merge_array;
   }

   
    /////////////////////////////////////////////////////////////
   //  Random double bwtween low and high
   /////////////////////////////////////////////////////////////
   public static double random_double_between (double low, double high)
    {
     double result=low+(high-low)*Math.random(); 
     return result;
    } 
	
	
   ////////////////////////////////////
   //  Log Gamma for Levy Flight Steps
   ////////////////////////////////////    
   public static double logGamma(double x)
   {
      double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
      double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
                       + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
                       +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
      return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
   }

   ////////////////////////////////////
   //  Gamma for Levy Flight Steps
   ////////////////////////////////////   
  public static double gamma(double x)
   {
     return Math.exp(logGamma(x));
   }
 
   ////////////////////////////////////
   //  Levy Flight Steps
   ////////////////////////////////////  
   public static double LevyFlight()
    {
      double beta=3/2;
      double sigma=Math.pow(gamma(1+beta)*Math.sin(Math.PI*beta/2)/(gamma((1+beta)/2)*beta*Math.pow(2,((beta-1)/2))),(1/beta));
      double u,v;
      Random randn = new Random();

       u=randn.nextGaussian()*sigma;
       v=randn.nextGaussian();

       //System.out.println (u);
       //System.out.println (v);
       double step= u/(Math.pow(Math.abs(v),(1/beta))) ;

      return Math.round(step);
    }

 
   /////////////////////////////////////////////////////////////
   //  Display any string array list for debugging purposes
   /////////////////////////////////////////////////////////////
   public final static void display_list (String title, ArrayList<String> list)
    {
     int i=0;
     System.out.println (title);
     for (Iterator it = list.iterator(); it.hasNext(); )
      {
        String s = (String)it.next();  // Downcasting is required pre Java 5.
        System.out.println (s);
        i++;
      }
    }
	
	/////////////////////////////////////////////////////////////
    //    Calculate standard deviation
    /////////////////////////////////////////////////////////////
    public static double calculate_stdv(double[] array_of_interests, double mean)
    {
        double sum=0.0;
        double variance=0.0;
        for(int i=0;i<array_of_interests.length;i++)
        {
            sum += Math.pow((array_of_interests[i]-mean), 2);
        }
        variance = sum/array_of_interests.length;
        return Math.sqrt(variance);
    }
	
    public static double roundTwoDecimals(double d) 
    {
        DecimalFormat two_d_form = new DecimalFormat("#.##");	
        return Double.valueOf(two_d_form.format(d));
    }
	
 }
