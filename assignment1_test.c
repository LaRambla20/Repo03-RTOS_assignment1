//compile with: g++ -pthread assignment1.c -o assignment1
//run with: sudo ./assignment1

// The first assignment consists in scheduling four tasks, each of which interacts with one or more 
// shared global variable. The access to this variables is protected by three semaphores handled with
// the priority ceiling protocol.
// This script shows one possible solution to this problem.

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <cstring>


//code that periodic tasks executes when implemented by threads
void task1_code( );
void task2_code( );
void task3_code( );
void task4_code( );

//characteristic function of the threads, only for timing and synchronization 
//between periodic tasks
void *task1( void *);
void *task2( void *);
void *task3( void *);
void *task4( void *);

// function to waste time (more specifically to increase the tasks computational time)
void waste_time( );

// initialization of mutexes for protecting the access shared global variables between tasks
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex3 = PTHREAD_MUTEX_INITIALIZER;

#define INNERLOOP 500 // global variables used in the function: waste_time()
#define OUTERLOOP 120

#define NPERIODICTASKS 4 // number of periodic tasks
#define NTASKS NPERIODICTASKS // total number of tasks
#define NTEST 30 // number of tests to evaluate the worst case execution time and the maximum blocking time 

// shared global variables:
int T1T2; // written by task 1, read by task 2
int T1T4; // written by task 1, read by task 4
int T2T3; // written by task 2, read by task 3


long int periods[NTASKS]; // array that stores as long ints the periods of the tasks
struct timespec next_arrival_time[NTASKS]; // array that stores as timespec structs the next arrival time of the tasks
double WCET[NTASKS]; // array that stores as doubles the Worst Case Execution Time of the tasks
double ET_current; // variable of type double that stores the current Execution Time of the considered task
double MBT[NTASKS]; // array that stores as doubles the Maximum Blocking Time for the tasks
double BT_z21_current; // variable that stores the Blocking Time related to the z21 critical section (first critical section of the second task)
double BT_z31_current; // variable that stores the Blocking Time related to the z31 critical section (first critical section of the third task)
double BT_z41_current; // variable that stores the Blocking Time related to the z41 critical section (first critical section of the fourth task)
double U[NTASKS]; // array that stores as doubles the Utilization Factors of the tasks
double Ulub[NTASKS]; // array that stores as doubles the lower-upper-bound Utilization Factors of the tasks
pthread_attr_t attributes[NTASKS]; // array that stores as pthread_attr_t the attributes related to each task (attributes contain information about the scheduling policy adopted)
pthread_mutexattr_t mutexattributes; //  variable of type pthread_mutexattr_t that stores the attributes related to the mutexes (mutexes attributes contain information about the protocol adopted to handle them)
pthread_t thread_id[NTASKS]; // array that stores as pthread_t the ids of the tasks
struct sched_param parameters[NTASKS]; // array that stores as sched_param structs the parameters related to each task (parameters conntain information about the priority of the tasks)
int missed_deadlines[NTASKS]; // array that stores as ints the number of missed deadlines for the tasks



int main()
{
  	// set task periods in nanoseconds
	// the first task has period 80 millisecond
	// the second task has period 100 millisecond
	// the third task has period 160 millisecond
    // the fourth task has period 200 millisecond
	// they are ordered according to their priority in RM scheduling algorithm (the priority is inversely proportional to the period: the higher the period, the lower the priority); 
	periods[0]= 80000000; //period of the first task in nanoseconds
  	periods[1]= 100000000; //in nanoseconds
  	periods[2]= 160000000; //in nanoseconds
    periods[3]= 200000000; //in nanoseconds

	//this is not strictly necessary, but it is convenient to
	//assign a name to the maximum and the minimum priority in the
	//system. We call them priomin and priomax.

  	struct sched_param priomax; // variable to store the maximum possible priority using the sched_FIFO
  	priomax.sched_priority=sched_get_priority_max(SCHED_FIFO); // function to get the maximum possible priority using the sched_FIFO
  	struct sched_param priomin; // variable to store the maximum possible priority using the sched_FIFO
  	priomin.sched_priority=sched_get_priority_min(SCHED_FIFO); // function to get the minimum possible priority using the sched_FIFO

	// Check that the main thread is executed with superuser privileges and 
	// then set the maximum priority to the current/main thread (this is because it is required to be
  	// superuser in order to do that)

	if (getuid() == 0)
    	pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax);


	// since blocking of tasks may occur because of the presence of mutexes, for every task (except for the one with the least priority), in addition to the computational time,
	// the maximum blocking time due to tasks with lower priority should be computed. Using the priority ceiling protocol, every task can be blocked by a lower priority task 
	// for at maximum 1 critical section. So the maximum blocking time is nothing more than the longest duration among the durations of critical sections that can block the 
	// task at issue. As far as the task with the least priority is concerned, its maximum blocking time is 0 beacuse there is no lower priority task that can block it

  	// execute all tasks in standalone modality NTEST times in order to measure the execution times and the maximum blocking times.
	// Based on the results, update the Worst Case Execution Time and the Maximum Blocking Time for each tasks.
	// The computed WCETs and MBTs are just estimates of the real values, since, to obtain them, each task should be executed on its own
	// for many more times than NTESTS.

    printf("TESTING\n\n");
    fflush(stdout);
    
    memset(WCET, 0, NTASKS);
	memset(MBT, 0, NTASKS);
 	int i;
    int j;
    for (i = 0; i < NTEST; i++){
        for (j =0; j < NTASKS; j++)
        {

            struct timespec time_1, time_2;
            clock_gettime(CLOCK_REALTIME, &time_1); // store the current time of the day in the time_1 variable before the j-th tasks starts computing
            // the current time is retrieved in seconds and nanoseconds

            if (j==0)
                task1_code();
            if (j==1)
                task2_code();
            if (j==2)
                task3_code();
            if (j==3)
                task4_code();

            clock_gettime(CLOCK_REALTIME, &time_2); // store the current time of the day in the time_2 variable after the j-th tasks finished computing

            ET_current= 1000000000*(time_2.tv_sec - time_1.tv_sec) + (time_2.tv_nsec-time_1.tv_nsec); // compute the current Execution Time of the considered task
            if (ET_current > WCET[j]){ // update the Worst Case Execution Time of the considered task if the current Execution Time is higher than the previous WCET stored
                WCET[j] = ET_current;
                // if(j==1){
                //     printf("the new value is: %lf\n", WCET[j]);
                // }
            }

        }

		// Retrieve the current Blocking Times of the critical sections related to each task and:

		//Update task1 Maximum Blocking Time if the value of the highest current Blocking Time (related to task1) is greater than the previous MBT stored
		if (BT_z21_current > BT_z41_current){
			if(BT_z21_current > MBT[0]){
				MBT[0] = BT_z21_current;
			}
		}
		else{
			if(BT_z41_current > MBT[0]){
				MBT[0] = BT_z41_current;
			}
		}

		//Update task2 Maximum Blocking Time if the value of the highest current Blocking Time (related to task2) is greater than the previous MBT stored
		if (BT_z31_current > BT_z41_current){
			if(BT_z31_current > MBT[1]){
				MBT[1] = BT_z31_current;
			}
		}
		else{
			if(BT_z41_current > MBT[1]){
				MBT[1] = BT_z41_current;
			}
		}

		//Update task3 Maximum Blocking Time if the value of the highest current Blocking Time (related to task3) is greater than the previous MBT stored
		if(BT_z41_current > MBT[2]){
			MBT[2] = BT_z41_current;
		}

    }


    for(i = 0; i < NTASKS; i++){
        printf("\nWorst Case Execution Time [%d]=%f \n", i, WCET[i]);
        fflush(stdout);
		printf("\nMaximum Blocking Time [%d]=%f \n", i, MBT[i]);
		fflush(stdout);
    }


	memset(U, 0, NTASKS);

	for (i =0; i < NTASKS; i++)
    {
		// compute U[i] (Utilization Factor of the i-th task) by considering also its Maximum Blocking Time
		for (j = 0; j <= i; j++){
			U[i] = U[i] + (WCET[j]/periods[j]);
		}
		U[i] = U[i] + MBT[i]/periods[i];

		// compute Ulub up to the i-th task considered
		Ulub[i] = (i+1)*(pow(2.0,(1.0/(i+1))) -1); // i+1 because, since i is an array index, it starts from 0, so the correspondent task number is i+1

		//check the sufficient conditions: if they are not satisfied, exit because it means that nothing can be said about the schedulability of the tasks with RM
		if (U[i] > Ulub[i])
    	{
      	printf("\n U[%d]=%lf Ulub[%d]=%lf Non schedulable Task Set", i, U[i], i, Ulub[i]);
		fflush(stdout);
      	return(-1);
    	}
  		printf("\n U[%d]=%lf Ulub[%d]=%lf Scheduable Task Set", i, U[i], i, Ulub[i]);
  		fflush(stdout);
		
    }

	printf("\n\n");
	fflush(stdout);
	
  	sleep(5);

  	// set the minimum priority to the current thread: this is now required because 
	// it will be assigned higher priority to periodic threads that will be soon created
  	if (getuid() == 0)
    	pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);

  
  	// set the attributes and the parameters for each task in order to specify the scheduling algorithm and the priorities
  	for (i =0; i < NPERIODICTASKS; i++)
    {
		//initialize the attribute structure of the i-th task
		pthread_attr_init(&(attributes[i])); 

		//set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		//not inherited from the main thread (pthread_attr_setinheritsched) 
		pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
		pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO); 

		//properly set the parameters to assign the priority inversely proportional to the period
		parameters[i].sched_priority = priomin.sched_priority+NTASKS - i;

		//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
		pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    	
	}


	// set the attributes of the semaphores in order to specify that they should be handeled with the priority ceiling protocol

	//initialize the attributes structure of the mutexes
	pthread_mutexattr_init(&mutexattributes);
		
	// set the attributes to set the PTHREAD_PRIO_PROTECT protocol (pthread_mutexattr_setprotocol)
	pthread_mutexattr_setprotocol(&mutexattributes, PTHREAD_PRIO_PROTECT);

	// for each mutex: set the priority ceiling protocol, specifying the correspondent ceiling (maximum priority among the priorities of the tasks that 
	// use that semaphore) (pthread_mutexattr_setprioceiling), and initialize it

    // UNCOMMENT TO SOLVE THE DEADLOCK (here the priority ceiling protocol is set)
	// pthread_mutexattr_setprioceiling(&mutexattributes, parameters[0].sched_priority); // ceiling: priority of task1
	// pthread_mutex_init(&mutex1, &mutexattributes);
	// pthread_mutexattr_setprioceiling(&mutexattributes, parameters[0].sched_priority); // ceiling: priority of task1
	// pthread_mutex_init(&mutex2, &mutexattributes);
	// pthread_mutexattr_setprioceiling(&mutexattributes, parameters[1].sched_priority); // ceiling: priority of task2
	// pthread_mutex_init(&mutex3, &mutexattributes);
    // UNCOMMENT TO SOLVE THE DEADLOCK


	//declare the array to contain the return values of pthread_create	
  	int iret[NTASKS];


 	// evaluate the second arrival time for each task, that is the end of the first period and the beginning of the second one
	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1); // store in the time_1 variable the initial time (before all the threads are run)
 
  	for (i = 0; i < NPERIODICTASKS; i++)
    {
        // first the current time in nanoseconds plus the period (that is in nanoseconds) is stored in a variable: in this way it can be obtained a result higher than 
        // 1000000000 nsec (that is: higher than 1 sec -> not higher than 2 or more seconds because those more seconds are contained in the "seconds field" of the struct).
        // So, in order to evaluate the residual nanoseconds till next_arrival_time, the additional number of nanoseconds with respect to 1000000000 must be obtained 
		// -> that's why the module is evaluated
        long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i]; 

        next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
        next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000; //convert nanoseconds into seconds
        missed_deadlines[i] = 0;
	}
	
	// create all threads(pthread_create)
	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
    iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);

  	// join all threads (pthread_join) / or rather make the main thread wait for the completion of the 4 threads
	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
    pthread_join( thread_id[3], NULL);


  	// print the number of missing deadlines for each task
  	for (i = 0; i < NTASKS; i++)
    {
    	printf ("\nMissed Deadlines Task %d=%d", i, missed_deadlines[i]); 
		fflush(stdout);
    }

	// destroy the attributes of the semaphores previously defined before exiting
	pthread_mutexattr_destroy(&mutexattributes);

  	exit(0);
}


// FUNCTIONS


// code executed by task1 (application dependent (e.g. acquiring pictures))
void task1_code()
{
	//print the id of the current task
  	printf("\x1b[31m" " 1[ " "\x1b[0m"); fflush(stdout);

	//call the function waste_time() 2 times to increase the computational time of the task
	int k;
	for (k = 0; k < 2; k++){
	waste_time();
	}


	pthread_mutex_lock(&mutex1); //lock the semaphore mutex1
	T1T2 = 1; // write on the T1T2 shared global variable
	pthread_mutex_unlock(&mutex1); //unlock the semaphore mutex1

	printf("(wT1T2:%d)", 1);
	fflush(stdout);

    pthread_mutex_lock(&mutex3);
    waste_time();
	pthread_mutex_lock(&mutex2); //lock the semaphore mutex2
	T1T4 = 1; // write on the T1T4 shared global variable
	pthread_mutex_unlock(&mutex2); //unlock the semaphore mutex2
    pthread_mutex_unlock(&mutex3);

	printf("(wT1T4:%d)", 1);
	fflush(stdout);
  
  	//print the id of the current task
  	printf("\x1b[31m" " ]1 " "\x1b[0m"); fflush(stdout);
}

//thread code for task1 (used for temporization/timing and synchronization between tasks, and to specify that all threads should run on the same core/processor/CPU)
void *task1( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset); // make the first thread run on the first CPU

   	//execute the task at maximum one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	for (i=0; i < 100; i++)
    {
      	//execute the task1 code
		task1_code();

		//check if the task missed a deadline during the i-th execution
        struct timespec time_t1;
	    clock_gettime(CLOCK_REALTIME, &time_t1); // store in the time_t1 variable the time at the end of the i-th execution of the task

		// if the current time is greater than the next arrival time previously evaluated, increment the number of deadlines missed by the task
        if((time_t1.tv_sec>next_arrival_time[0].tv_sec)||((time_t1.tv_sec==next_arrival_time[0].tv_sec)&&(time_t1.tv_nsec>next_arrival_time[0].tv_nsec))){
            missed_deadlines[0] = missed_deadlines[0]+1;
        }


		// sleep until a known given time (that is the beginning of the next period: next_arrival_time[0] in this case)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL); 

		// compute the end of the new period for the next iteration (the first iteration of the loop is carried out thanks to the fact that, before the creation of the threads, 
		// the beginning of the second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}



// code executed by task2
void task2_code()
{
	//print the id of the current task
  	printf("\x1b[33m" " 2[ " "\x1b[0m"); fflush(stdout);

	//call the function waste_time() 3 times to increase the computational time of the task
	int k;
    for (k = 0; k < 3; k++){
	waste_time();
	}

	// evaluate the Blocking Time related to the critical section z21
	struct timespec time_z21_beg;
	struct timespec time_z21_end;
	int read_value;
	clock_gettime(CLOCK_REALTIME, &time_z21_beg); // store in the time_z21_beg variable the time at the beginning of the critical section z21
	pthread_mutex_lock(&mutex1); //lock the semaphore mutex1
	read_value = T1T2; // read from T1T2 shared global variable
	pthread_mutex_unlock(&mutex1); //unlock the semaphore mutex1
	clock_gettime(CLOCK_REALTIME, &time_z21_end); // store in the time_z21_end variable the time at the end of the critical section z21

	printf("(rT1T2:%d)", read_value);
	fflush(stdout);

	// store the Blocking Time related to the critical section z21 in the global variable BT_z21_current
	BT_z21_current= 1000000000*(time_z21_end.tv_sec - time_z21_beg.tv_sec) + (time_z21_end.tv_nsec-time_z21_beg.tv_nsec);

    pthread_mutex_lock(&mutex2);
    waste_time();
    waste_time();
	pthread_mutex_lock(&mutex3); //lock the semaphore mutex3
	T2T3 = 2; // write on the T2T3 shared global variable
	pthread_mutex_unlock(&mutex3); //lock the semaphore mutex3
    pthread_mutex_unlock(&mutex2); 

	printf("(wT2T3:%d)", 2);
	fflush(stdout);

	//print the id of the current task
  	printf("\x1b[33m" " ]2 " "\x1b[0m"); fflush(stdout);
}

//thread code for task2 (used for temporization/timing and synchronization between tasks, and to specify that all threads should run on the same core/processor/CPU)
void *task2( void *ptr )
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	//execute the task at maximum one hundred times... it should be an infinite loop (too dangerous)
	int i=0;
  	for (i=0; i < 100; i++)
    {
		
        //execute the task2 code
		task2_code();

        //check if the task missed a deadline during the i-th execution
        struct timespec time_t2; 
	    clock_gettime(CLOCK_REALTIME, &time_t2); // store in the time_t2 variable the time at the end of the i-th execution of the task

		// if the current time is greater than the next arrival time previously evaluated, increment the number of deadlines missed by the task
        if((time_t2.tv_sec>next_arrival_time[1].tv_sec)||((time_t2.tv_sec==next_arrival_time[1].tv_sec)&&(time_t2.tv_nsec>next_arrival_time[1].tv_nsec))){
              missed_deadlines[1] = missed_deadlines[1]+1;
        }

        // sleep until a known given time (that is the beginning of the next period: next_arrival_time[1] in this case)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);


        // compute the end of the new period for the next iteration (the first iteration of the loop is carried out thanks to the fact that, before the creation of the threads, 
		// the beginning of the second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1]; 
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;

    	}
}



// code executed by task3
void task3_code()
{
	//print the id of the current task
  	printf("\x1b[32m" " 3[ " "\x1b[0m"); fflush(stdout);

	//call the function waste_time() 5 times to increase the computational time of the task
	int k;
    for (k = 0; k < 5; k++){
	waste_time();
	}

	// evaluate the Blocking Time related to the critical section z31
	struct timespec time_z31_beg;
	struct timespec time_z31_end;
	int read_value;
	clock_gettime(CLOCK_REALTIME, &time_z31_beg); // store in the time_z31_beg variable the time at the beginning of the critical section z31
	pthread_mutex_lock(&mutex3); //lock the semaphore mutex3
	read_value = T2T3; // read from the T2T3 shared global variable
	pthread_mutex_unlock(&mutex3); //unlock the semaphore mutex3
	clock_gettime(CLOCK_REALTIME, &time_z31_end); // store in the time_z31_end variable the time at the end of the critical section z31

	printf("(rT2T3:%d)", read_value);
	fflush(stdout);

	// store the Blocking Time related to the critical section z31 in the global variable BT_z31_current
	BT_z31_current= 1000000000*(time_z31_end.tv_sec - time_z31_beg.tv_sec) + (time_z31_end.tv_nsec-time_z31_beg.tv_nsec);

	//print the id of the current task
  	printf("\x1b[32m" " ]3 " "\x1b[0m"); fflush(stdout);
}

//thread code for task3 (used for temporization/timing and synchronization between tasks, and to specify that all threads should run on the same core/processor/CPU)
void *task3( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	//execute the task at maximum one hundred times... it should be an infinite loop (too dangerous)
	int i=0;
  	for (i=0; i < 100; i++)
    {

        //execute the task3 code
		task3_code();

        //check if the task missed a deadline during the i-th execution
        struct timespec time_t3; 
	    clock_gettime(CLOCK_REALTIME, &time_t3); // store in the time_t3 variable the time at the end of the i-th execution of the task

		// if the current time is greater than the next arrival time previously evaluated, increment the number of deadlines missed by the task
        if((time_t3.tv_sec>next_arrival_time[2].tv_sec)||((time_t3.tv_sec==next_arrival_time[2].tv_sec)&&(time_t3.tv_nsec>next_arrival_time[2].tv_nsec))){
              missed_deadlines[2] = missed_deadlines[2]+1;
        }

        // sleep until a known given time (that is the beginning of the next period: next_arrival_time[2] in this case)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);

		// compute the end of the new period for the next iteration (the first iteration of the loop is carried out thanks to the fact that, before the creation of the threads, 
		// the beginning of the second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
		next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
	
	}
}



// code executed by task4
void task4_code()
{
	//print the id of the current task
  	printf("\x1b[34m" " 4[ " "\x1b[0m"); fflush(stdout);

	//call the function waste_time() 6 times to increase the computational time of the task
	int k;
    for (k = 0; k < 6; k++){
	waste_time();
	}

	// evaluate the Blocking Time related to the critical section z41
	struct timespec time_z41_beg;
	struct timespec time_z41_end; 
	int read_value;
	clock_gettime(CLOCK_REALTIME, &time_z41_beg); // store in the time_z41_beg variable the time at the beginning of the critical section z41
	pthread_mutex_lock(&mutex2); //lock the semaphore mutex2
	read_value = T1T4; //write on the T1T4 shared global variable
	pthread_mutex_unlock(&mutex2); //unlock the semaphore mutex1
	clock_gettime(CLOCK_REALTIME, &time_z41_end); // store in the time_z41_end variable the time at the end of the critical section z41

	printf("(rT1T4:%d)", read_value);
	fflush(stdout);

	// store the Blocking Time related to the critical section z21 in the global variable BT_z41_current
	BT_z41_current= 1000000000*(time_z41_end.tv_sec - time_z41_beg.tv_sec) + (time_z41_end.tv_nsec-time_z41_beg.tv_nsec);

	//print the id of the current task
  	printf("\x1b[34m" " ]4 " "\x1b[0m"); fflush(stdout);
}

//thread code for task4 (used for temporization/timing and synchronization between tasks, and to specify that all threads should run on the same core/processor/CPU)
void *task4( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	//execute the task at maximum one hundred times... it should be an infinite loop (too dangerous)
	int i=0;
  	for (i=0; i < 100; i++)
    {

        //execute the task4 code
		task4_code();

        //check if the task missed a deadline during the i-th execution
        struct timespec time_t4;
	    clock_gettime(CLOCK_REALTIME, &time_t4); // store in the time_t4 variable the time at the end of the i-th execution of the task

		// if the current time is greater than the next arrival time previously evaluated, increment the number of deadlines missed by the task
        if((time_t4.tv_sec>next_arrival_time[3].tv_sec)||((time_t4.tv_sec==next_arrival_time[3].tv_sec)&&(time_t4.tv_nsec>next_arrival_time[3].tv_nsec))){
              missed_deadlines[3] = missed_deadlines[3]+1;
        }

        // sleep until a known given time (that is the beginning of the next period: next_arrival_time[3] in this case)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);

		// compute the end of the new period for the next iteration (the first iteration of the loop is carried out thanks to the fact that, before the creation of the threads, 
		// the beginning of the second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
		next_arrival_time[3].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec= next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
	
	}
}




// function that implements a loop innested into another one and the aim of which is to waste a certain amount of time
void waste_time(){
	int i,j;
	double uno;
  	for (i = 0; i < OUTERLOOP; i++) // OUTERLOOP is a global constant
    {
      	for (j = 0; j < INNERLOOP; j++) // INNERLOOP is a global constant
		{
			uno = rand()*rand()%10; // random operation
		}
    }
}