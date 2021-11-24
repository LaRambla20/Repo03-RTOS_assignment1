//compile with: g++ -pthread assignment1.c -o assignment1
//run with: sudo ./assignment1

//This exercise shows how to schedule threads with Rate Monotonic

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

// function to waste time (more specifically to increase the tasks' computational time)
void waste_time( );

// initialization of mutexes for protecting the reading/writing of shared global variables between tasks
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex3 = PTHREAD_MUTEX_INITIALIZER;

#define INNERLOOP 800 // global variables
#define OUTERLOOP 80

#define NPERIODICTASKS 4 // number of periodic tasks
#define NTASKS NPERIODICTASKS // total number of tasks

// shared global variables:
int T1T2; // written by task 1, read by task 2
int T1T4; // written by task 1, read by task 4
int T2T3; // written by task 2, read by task 3


long int periods[NTASKS]; // array of periods, stored as long ints, with NTASKS cells
struct timespec next_arrival_time[NTASKS]; // array of next_arrival_time, stored as structs timespec (struct with certain particular fields), with NTASKS cells
double WCET[NTASKS]; // array of Worst Case Execution Time, stored as doubles, with NTASKS cells
double MBT[NTASKS]; // array of the Maximum Blocking Time, stored as doubles, with NTASKS cells
double U[NTASKS]; // array of Utilization Factors, stored as doubles, with NTASKS cells
double Ulub[NTASKS]; // array of U lower upper bound, stored as doubles, with NTASKS cells
pthread_attr_t attributes[NTASKS]; // array of tasks attributes, stored as pthread_attr_t, with NTASKS cells
pthread_mutexattr_t mutexattributes; //  variable of type pthread_mutexattr_t that stores the attributes of the mutexes
pthread_t thread_id[NTASKS]; // array of thread_id, stored as pthread_t, with NTASKS cells
struct sched_param parameters[NTASKS]; // array of parameters of the scheduling algorithm, stored as structs sched_param, with NTASKS cells
int missed_deadlines[NTASKS]; // array of number of missed deadines for each task, stored as ints, with NTASKS cells



int main()
{
  	// set task periods in nanoseconds
	//the first task has period 80 millisecond
	//the second task has period 100 millisecond
	//the third task has period 160 millisecond
    //the fourth task has period 200 millisecond
	//they are ordered according to their priority in RM scheduling algorithm (the priority is inversely proportional to the period: the higher the period, the lower the priority); 
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
	// then set the maximum priority to the current/main thread (this is because you are required to be
  	// superuser in order to do that)

	if (getuid() == 0)
    	pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax);


  	// execute all tasks in standalone modality in order to measure execution times
  	// (using clock_gettime). Use the computed values to update the worst case execution
  	// time of each task.

 	int i;
  	for (i =0; i < NTASKS; i++)
    {
		// initialize time_1 and time_2 required to read the clock
		struct timespec time_1, time_2; // two variables stored as structs timespec (struct with certain particular fields)
		clock_gettime(CLOCK_REALTIME, &time_1); // store the current time of the day in the time_1 variable before the i-th tasks starts computing
		// the current time is retrieved in seconds and nanoseconds

		// we should execute each task on their own for many times to compute the real WCET
		// here we do that just once, so it is not an accurate estimate

		if (i==0)
		    task1_code();
      	if (i==1)
		    task2_code();
      	if (i==2)
		    task3_code();
        if (i==3)
		    task4_code();

		clock_gettime(CLOCK_REALTIME, &time_2); // store the current time of the day in the time_2 variable after the i-th tasks finished computing
		// the current time is retrieved in seconds and nanoseconds

		// compute an estimate of the Worst Case Execution Time (because we executed each task just once)

		WCET[i]= 1000000000*(time_2.tv_sec - time_1.tv_sec) + (time_2.tv_nsec-time_1.tv_nsec); // seconds must be converted in nanoseconds: we work in nanoseconds
      	printf("\nWorst Case Execution Time [%d]=%f \n", i, WCET[i]);
		fflush(stdout);

    } // end of the for

	// since blocking of tasks may occur because of the presence of mutexes, for every task (except for the one with the least priority) the maximum blocking time
	// due to tasks with lower priority should be computed. Using the priority ceiling protocol, every task can be blocked by a lower priority task for at maximum
	// 1 critical section. So the maximum blocking time is nothing more than the longest duration among the durations of critical sections that can block the task at issue.
	// It can be proven that each task (except for the one with the least priority) can only be blocked by tasks with lower priority when these are reading from the 
	// shared variables, so the maximum blocking time corresponds to the time taken to read. As far as the task with the least priority is concerned, its maximum blocking
	// time is 0 beacuse there is no lower priority task that can block it

	int read_value; // variable to read the shared global variables

	// initialize time_1 and time_2 required to read the clock
	struct timespec time_1, time_2;
	clock_gettime(CLOCK_REALTIME, &time_1);

	read_value = T1T2;

	clock_gettime(CLOCK_REALTIME, &time_2);

	for (i = 0; i < NTASKS; i++){
		if (i != NTASKS-1){
			MBT[i] = 1000000000*(time_2.tv_sec - time_1.tv_sec) + (time_2.tv_nsec-time_1.tv_nsec); // seconds must be converted in nanoseconds: we work in nanoseconds
		}
		else{
			MBT[i] = 0;
		}
    	printf("\nMaximum Blocking Time [%d]=%f \n", i, MBT[i]);
		fflush(stdout);
	}



	memset(U, 0, NTASKS); // initialize all th elements of the vector to 0

	int j;
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
	//we will assign higher priorities to periodic threads to be soon created
  	if (getuid() == 0)
    	pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);

  
  	// set the attributes of each task, including scheduling policy and priority
  	for (i =0; i < NPERIODICTASKS; i++)
    {
		//initialize the attribute structure of the i-th task
		pthread_attr_init(&(attributes[i])); // attributes contain information about the scheduling policy adopted

		//set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		//not inherited from the main thread (pthread_attr_setinheritsched) 
		pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
		pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO); // set the scheduling policy

		//properly set the parameters to assign the priority inversely proportional to the period
		parameters[i].sched_priority = priomin.sched_priority+NTASKS - i; // I take the minimum priority possible, 
        // I sum the total number of tasks and every time I subtract the number of the task at issue
		// because tasks are ordered in increasing period order 
        // (imagine that priority is an int. The higher the priority, the higher the int: if the minimum priority is 1, I add the number of tasks (3) and
		// each loop I substract the number of the task at issue (es task 2): the higher the the number of the task at issue, the lower the priority)

		// parameters contain information about tasks priority

		//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
		pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    	
	} // end of the for


	// set the attributes of the semaphores in order to specify that they should be handeled with the priority ceiling protocol

	//initialize the attributes structure of the mutexes
	pthread_mutexattr_init(&mutexattributes);
		
	// set the attributes to set the PTHREAD_PRIO_PROTECT protocol (pthread_mutexattr_setprotocol)
	pthread_mutexattr_setprotocol(&mutexattributes, PTHREAD_PRIO_PROTECT);

	// for each mutex: set the priority ceiling protocol, specifying the correspondent ceiling (maximum priority among the priorities of the tasks that 
	// use that semaphore) (pthread_mutexattr_setprioceiling), and initialize it
	pthread_mutexattr_setprioceiling(&mutexattributes, 1); // ceiling: 1 (priority of task1)
	pthread_mutex_init(&mutex1, &mutexattributes);
	pthread_mutexattr_setprioceiling(&mutexattributes, 1); // ceiling: 1 (priority of task1)
	pthread_mutex_init(&mutex2, &mutexattributes);
	pthread_mutexattr_setprioceiling(&mutexattributes, 2); // ceiling: 2 (priority of task2)
	pthread_mutex_init(&mutex3, &mutexattributes);


	//declare the array to contain the return values of pthread_create	
  	int iret[NTASKS];

	clock_gettime(CLOCK_REALTIME, &time_1); // store the current time of the day in the time_1 variable
	// the current time is retrieved in seconds and nanoseconds

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NPERIODICTASKS; i++)
    {
        // first we store the current time in nanoseconds plus the period (that is in nanoseconds) in a variable: this way I can obtain a result higher than 
        // 1000000000 nsec (that is: higher than 1 sec -> not higher than 2 or more seconds because that more seconds are contained in the "seconds field" of the struct).
        // So to evaluate the residual nanoseconds till next_arrival_time, I have to obtain the additional number of
        // nanoseconds with respect to 1000000000 -> that's why I evaluate the module
        long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i]; 

        //then we compute the end of the first period and beginning of the next one
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


// code executed by task 1 (application dependent e.g. acquiring pictures)
void task1_code()
{
	//print the id of the current task
  	printf(" 1[ "); fflush(stdout);

	int k;
	for (k = 0; k < 2; k++){
	//call the function waste_time() to increase the computational time of the task
	waste_time();
	}

	pthread_mutex_lock(&mutex1);
	T1T2 = 1;
	pthread_mutex_unlock(&mutex1);

	printf("(wT1T2:%d)", 1);
	fflush(stdout);

	pthread_mutex_lock(&mutex2);
	T1T4 = 1;
	pthread_mutex_unlock(&mutex2);

	printf("(wT1T4:%d)", 1);
	fflush(stdout);
  
  	//print the id of the current task
  	printf(" ]1 "); fflush(stdout);
}

//thread code for task_1 (used only for temporization / timing and synchronization between tasks)
// (and also to specify that all threads should run on the same core/processor/CPU)
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
      	//execute the application dependent code
		task1_code();

		//check if the task missed a deadline during the i-th execution

        struct timespec time_t1; // struct to store the current time in the characteristic function of the first thread: void *task1( void *), in order to check if some deadlines have been missed
	    clock_gettime(CLOCK_REALTIME, &time_t1); // store the current time of the day in the time_1 variable

        if((time_t1.tv_sec>next_arrival_time[0].tv_sec)||((time_t1.tv_sec==next_arrival_time[0].tv_sec)&&(time_t1.tv_nsec>next_arrival_time[0].tv_nsec))){
              missed_deadlines[0] = missed_deadlines[0]+1;
        }


		// sleep until the end of the current period (which is also the start of the new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL); // make the thread sleep until a known given time (that is the beginning
		// of the next period: next_arrival_time[0])

		// compute the end of the current period for the next iteration

        //just before exiting the loop, I update the variable next_arrival_time[0] (the first iteration of the loop is carried out thanks to the fact that, before the loop, the beginning of the
        // second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0]; // each time the thread is executed the subsequent arrival time (that is the beginning of the subsequent period)
		// is evaluated as seen before (N.B. before only the beginning of the second period was evaluated)
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}



void task2_code()
{
	//print the id of the current task
  	printf(" 2[ "); fflush(stdout);

	int k;
    for (k = 0; k < 3; k++){
	//call the function waste_time() to increase the computational time of the task
	waste_time();
	}

	pthread_mutex_lock(&mutex3);
	T2T3 = 2;
	pthread_mutex_unlock(&mutex3);

	printf("(wT2T3:%d)", 2);
	fflush(stdout);


	int read_value;
	pthread_mutex_lock(&mutex1);
	read_value = T1T2;
	pthread_mutex_unlock(&mutex1);

	printf("(rT1T2:%d)", read_value);
	fflush(stdout);

	//print the id of the current task
  	printf(" ]2 "); fflush(stdout);
}


void *task2( void *ptr )
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
  	for (i=0; i < 100; i++)
    {
		
        // execute the application dependent code
		task2_code();

        //check if the task missed a deadline during the i-th execution

        struct timespec time_t2; // struct to store the current time in the characteristic function of the first thread: void *task1( void *), in order to check if some deadlines have been missed
	    clock_gettime(CLOCK_REALTIME, &time_t2); // store the current time of the day in the time_1 variable

        if((time_t2.tv_sec>next_arrival_time[1].tv_sec)||((time_t2.tv_sec==next_arrival_time[1].tv_sec)&&(time_t2.tv_nsec>next_arrival_time[1].tv_nsec))){
              missed_deadlines[1] = missed_deadlines[1]+1;
        }

        // sleep until the end of the current period (which is also the start of the new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL); // make the thread sleep until a known given time (that is the beginning
		// of the next period: next_arrival_time[0])


        // compute the end of the current period for the next iteration

        //just before exiting the loop, I update the variable next_arrival_time[0] (the first iteration of the loop is carried out thanks to the fact that, before the loop, the beginning of the
        // second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1]; // each time the thread is executed the subsequent arrival time (that is the beginning of the subsequent period)
		// is evaluated as seen before (N.B. before only the beginning of the second period was evaluated)
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;

    	}
}



void task3_code()
{
	//print the id of the current task
  	printf(" 3[ "); fflush(stdout);

	int k;
    for (k = 0; k < 5; k++){
	//call the function waste_time() to increase the computational time of the task
	waste_time();
	}

	int read_value;
	pthread_mutex_lock(&mutex3);
	read_value = T2T3;
	pthread_mutex_unlock(&mutex3);

	printf("(rT2T3:%d)", read_value);
	fflush(stdout);

	//print the id of the current task
  	printf(" ]3 "); fflush(stdout);
}

void *task3( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
  	for (i=0; i < 100; i++)
    {

        //execute the application dependent code
		task3_code();

        //check if the task missed a deadline during the i-th execution

        struct timespec time_t3; // struct to store the current time in the characteristic function of the first thread: void *task1( void *), in order to check if some deadlines have been missed
	    clock_gettime(CLOCK_REALTIME, &time_t3); // store the current time of the day in the time_1 variable

        if((time_t3.tv_sec>next_arrival_time[2].tv_sec)||((time_t3.tv_sec==next_arrival_time[2].tv_sec)&&(time_t3.tv_nsec>next_arrival_time[2].tv_nsec))){
              missed_deadlines[2] = missed_deadlines[2]+1;
        }

        // sleep until the end of the current period (which is also the start of the new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);// make the thread sleep until a known given time (that is the beginning
		// of the next period: next_arrival_time[0])

        // compute the end of the current period for the next iteration

        //just before exiting the loop, I update the variable next_arrival_time[0] (the first iteration of the loop is carried out thanks to the fact that, before the loop, the beginning of the
        // second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2]; // each time the thread is executed the subsequent arrival time (that is the beginning of the subsequent period)
		// is evaluated as seen before (N.B. before only the beginning of the second period was evaluated)
		next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
	
	}
}


void task4_code()
{
	//print the id of the current task
  	printf(" 4[ "); fflush(stdout);

	int k;
    for (k = 0; k < 6; k++){
	//call the function waste_time() to increase the computational time of the task
	waste_time();
	}

	int read_value;
	pthread_mutex_lock(&mutex2);
	read_value = T1T4;
	pthread_mutex_unlock(&mutex2);

	printf("(rT1T4:%d)", read_value);
	fflush(stdout);

	//print the id of the current task
  	printf(" ]4 "); fflush(stdout);
}

void *task4( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
  	for (i=0; i < 100; i++)
    {

        //execute the application dependent code
		task4_code();

        //check if the task missed a deadline during the i-th execution

        struct timespec time_t4; // struct to store the current time in the characteristic function of the first thread: void *task1( void *), in order to check if some deadlines have been missed
	    clock_gettime(CLOCK_REALTIME, &time_t4); // store the current time of the day in the time_1 variable

        if((time_t4.tv_sec>next_arrival_time[3].tv_sec)||((time_t4.tv_sec==next_arrival_time[3].tv_sec)&&(time_t4.tv_nsec>next_arrival_time[3].tv_nsec))){
              missed_deadlines[3] = missed_deadlines[3]+1;
        }

        // sleep until the end of the current period (which is also the start of the new one)
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);// make the thread sleep until a known given time (that is the beginning
		// of the next period: next_arrival_time[0])

        // compute the end of the current period for the next iteration

        //just before exiting the loop, I update the variable next_arrival_time[0] (the first iteration of the loop is carried out thanks to the fact that, before the loop, the beginning of the
        // second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3]; // each time the thread is executed the subsequent arrival time (that is the beginning of the subsequent period)
		// is evaluated as seen before (N.B. before only the beginning of the second period was evaluated)
		next_arrival_time[3].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec= next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
	
	}
}




void waste_time(){
	int i,j;
	double uno;
  	for (i = 0; i < OUTERLOOP; i++)
    {
      	for (j = 0; j < INNERLOOP; j++)
		{
			uno = rand()*rand()%10;
		}
    }
}
