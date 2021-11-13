//compile with: g++ -lpthread <sourcename> -o <executablename>

//This exercise shows how to schedule threads with Rate Monotonic

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>

//code that periodic tasks executes when implemented by threads
void task1_code( );
void task2_code( );
void task3_code( );

//code of aperiodic tasks (if any)
void task4_code( );
void task5_code( );
void task6_code( );

//characteristic function of the threads, only for timing and synchronization 
//between periodic tasks
void *task1( void *);
void *task2( void *);
void *task3( void *);

//aperiodic tasks (if any)
void *task4( void *);
void *task5( void *);
void *task6( void *);

// initialization of mutexes and conditions (only for aperiodic scheduling)
pthread_mutex_t mutex_task_4 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_task_5 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_task_6 = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond_task_4 = PTHREAD_COND_INITIALIZER;
pthread_cond_t cond_task_5 = PTHREAD_COND_INITIALIZER;
pthread_cond_t cond_task_6 = PTHREAD_COND_INITIALIZER;

#define INNERLOOP 1000 // global variables
#define OUTERLOOP 800

#define NPERIODICTASKS 3 // number of periodic tasks
#define NAPERIODICTASKS 3 // number of aperiodic tasks
#define NTASKS NPERIODICTASKS + NAPERIODICTASKS // total number of tasks

long int periods[NTASKS]; // array of periods, stored as long ints, with NTASKS cells
struct timespec next_arrival_time[NTASKS]; // array of next_arrival_time, stored as structs timespec (struct with certain particular fields), with NTASKS cells
double WCET[NTASKS]; // array of Worst Case Execution Time, stored as dubles, with NTASKS cells
pthread_attr_t attributes[NTASKS]; // array of attributes, stored as pthread_attr_t, with NTASKS cells
pthread_t thread_id[NTASKS]; // array of thread_id, stored as pthread_t, with NTASKS cells
struct sched_param parameters[NTASKS]; // array of parameters of the scheduling algorithm, stored as structs sched_param, with NTASKS cells
int missed_deadlines[NTASKS]; // array of number of missed deadines for each task, stored as ints, with NTASKS cells

int
main()
{
  	// set task periods in nanoseconds
	//the first task has period 100 millisecond
	//the second task has period 200 millisecond
	//the third task has period 400 millisecond
	//you can already order them according to their priority; 
	//if not, you will need to sort them
	periods[0]= 100000000; //period of the first task in nanoseconds
  	periods[1]= 200000000; //in nanoseconds
  	periods[2]= 400000000; //in nanoseconds

  	//for aperiodic tasks we set the period equals to 0
  	periods[3]= 0; 
  	periods[4]= 0; 
    periods[5]= 0;

	//this is not strictly necessary, but it is convenient to
	//assign a name to the maximum and the minimum priority in the
	//system. We call them priomin and priomax.

  	struct sched_param priomax; // variable to store the maximum possible priority using the sched_FIFO
  	priomax.sched_priority=sched_get_priority_max(SCHED_FIFO); // function to get the maximum possible priority using the sched_FIFO
  	struct sched_param priomin; // variable to store the maximum possible priority using the sched_FIFO
  	priomin.sched_priority=sched_get_priority_min(SCHED_FIFO);

	// Check that the main thread is executed with superuser privileges and 
	// then set the maximum priority to the current/main thread (this is because you are required to be
  	// superuser in order to do that)

	if (getuid() == 0)
    	pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax);


  	// execute all tasks in standalone modality in order to measure execution times
  	// (use gettimeofday). Use the computed values to update the worst case execution
  	// time of each task.

 	int i;
  	for (i =0; i < NTASKS; i++)
    	{

		// initialize time_1 and time_2 required to read the clock
		struct timespec time_1, time_2; // two variables stored as structs timespec (struct with certain particular fields)
		clock_gettime(CLOCK_REALTIME, &time_1); // store the current time of the day in the time_1 variable before the i-th tasks starts computing
		// in seconds and nanoseconds

		// we should execute each task on their own for many times computing the real WCET
		// here we do that just once
		//periodic tasks

			if (i==0)
			task1_code();
      		if (i==1)
			task2_code();
      		if (i==2)
			task3_code();
      		
      		//aperiodic tasks
            if (i==3)
			task4_code();
      		if (i==4)
			task5_code();
            if (i==5)
			task6_code();


		clock_gettime(CLOCK_REALTIME, &time_2); // store the current time of the day in the time_2 variable after the i-th tasks finished computing
		// in seconds and nanoseconds

		// compute the Worst Case Execution Time 
		// (in a real case, we should repeat this many times under different conditions, in order to have reliable values
		// here we do that just once

			WCET[i]= 1000000000*(time_2.tv_sec - time_1.tv_sec) // seconds must be converted in nanoseconds: we work in nanoseconds
			       +(time_2.tv_nsec-time_1.tv_nsec);
      		printf("\nWorst Case Execution Time %d=%f \n", i, WCET[i]);

    	} // end of the for

    	// compute U
		double U = WCET[0]/periods[0]+WCET[1]/periods[1]+WCET[2]/periods[2]; // I could also have used a FOR loop

    	// compute Ulub by considering the fact that we have harmonic relationships between periods
		double Ulub = 1;
    	
	//if there are no harmonic relationships, use the following formula instead
	//double Ulub = NPERIODICTASKS*(pow(2.0,(1.0/NPERIODICTASKS)) -1);
	
	//check the sufficient conditions: if they are not satisfied, exit because it means that the task set is not schedulable with RM
	if (U > Ulub)
    	{
      		printf("\n U=%lf Ulub=%lf Non schedulable Task Set", U, Ulub);
      		return(-1);
    	}
  	printf("\n U=%lf Ulub=%lf Scheduable Task Set", U, Ulub);
  	fflush(stdout);
  	sleep(5);

  	// set the minimum priority to the current thread: this is now required because 
	//we will assign higher priorities to periodic threads to be soon created
	//pthread_setschedparam

  	if (getuid() == 0)
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);

  
  	// set the attributes of each task, including scheduling policy and priority
  	for (i =0; i < NPERIODICTASKS; i++)
    	{
		//initialize the attribute structure of task i
			pthread_attr_init(&(attributes[i])); // attributes contain information about the scheduling policy adopted

		//set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		//not inherited from the main thread (pthread_attr_setinheritsched) 
			pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
			pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO); // set the scheduling policy

		//properly set the parameters to assign the priority inversely proportional 
		//to the period
			parameters[i].sched_priority = priomin.sched_priority+NTASKS - i; // I take the minimum priority possible, I sum the total number of tasks and every time I subtract the number of the task at issue
			// because tasks are ordered in increasing period order (imagine that priority is an int. Higher the priority higher the int: if the minimum priority is 1, I add the number of tasks (3) and
			// each loop I substract the number of the task at issue (es 3): higher the the number of task at issue, lower the priority)

			// parameters contain information about tasks priority

		//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
			pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    	
		} // end of the for

 	// aperiodic tasks
    for (int i =NAPERIODICTASKS; i < NTASKS; i++)
    {
      	pthread_attr_init(&(attributes[i]));
      	pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO);

      	//set minimum priority (background scheduling)
      	parameters[i].sched_priority = 0;
      	pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    }

	//declare the variable to contain the return values of pthread_create	
  	int iret[NTASKS];

	//declare variables to read the current time
	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1); // store the current time of the day in the time_1 variable
	// in seconds and nanoseconds

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NPERIODICTASKS; i++)
    	{
		// first we store the current time in nanoseconds plus the period (that is in nanoseconds) in a variable: this way I can obtain a result higher than 
		// 1000000000 nsec (AKA higher than 1 sec -> not higher than 2 or more seconds because that more seconds are contained in the seconds field of the struct).
		// So to evaluate the residual nanoseconds till next_arrival_time, I have to obtain the additional number of
		// nanoseconds with respect to 1000000000 -> that's why I evaluate the module
		long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i]; 

		//then we compute the end of the first period and beginning of the next one
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000;
       		missed_deadlines[i] = 0;
    	}
	
	// create all threads(pthread_create)
	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
    iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);
  	iret[4] = pthread_create( &(thread_id[4]), &(attributes[4]), task5, NULL);
    iret[5] = pthread_create( &(thread_id[5]), &(attributes[5]), task6, NULL);

  	// join all threads (pthread_join) / AKA make the main thread wait for the completion of the 3 threads
	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);


  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NTASKS; i++)
    	{
      		printf ("\nMissed Deadlines Task %d=%d", i, missed_deadlines[i]); // print the number of missing deadlines for each task
			fflush(stdout);
    	}
  	exit(0);
}


// FUNCTIONS

//- PERIODIC TASKS


// code executed by task 1 (application dependent e.g. acquiring pictures)
void task1_code()
{
	//print the id of the current task
  	printf(" 1[ "); fflush(stdout);

	//this double loop with random computation is only required to waste time
	int i,j;
	double uno;
  	for (i = 0; i < OUTERLOOP; i++)
    	{
      		for (j = 0; j < INNERLOOP; j++)
		{
			uno = rand()*rand()%10;
    		}
  	}

  	// when the random variable uno=0, then aperiodic task 5 must
  	// be executed/AKA the signal to task 4 is sent
    if (uno == 0)
    {
      	printf(":ex(4)");fflush(stdout);
    // In theory, we should protect conditions using mutexes. However, in a real-time application, something undesirable may happen.
	// Indeed, when task1 takes the mutex and sends the condition, task4 is executed and is given the mutex by the kernel. Which means
	// that task1 (higher priority) would be blocked waiting for task4 to finish (lower priority). This is of course unacceptable,
	// as it would produced a priority inversion. For this reason, we are not putting mutexes here. A better solution should be found.

//      	pthread_mutex_lock(&mutex_task_4);
      	pthread_cond_signal(&cond_task_4);
//      	pthread_mutex_unlock(&mutex_task_4);
	}

  	// when the random variable uno=1, then aperiodic task 5 must
  	// be executed
    if (uno == 1)
    {
    	printf(":ex(5)");fflush(stdout);

	// See below why mutexes have been commented
//  		pthread_mutex_lock(&mutex_task_5);
    	pthread_cond_signal(&cond_task_5);
//  		pthread_mutex_unlock(&mutex_task_5);
    }

    if (uno == 2)
    {
    	printf(":ex(6)");fflush(stdout);

	// See below why mutexes have been commented
//  		pthread_mutex_lock(&mutex_task_6);
    	pthread_cond_signal(&cond_task_6);
//  		pthread_mutex_unlock(&mutex_task_6);
    }
 

  
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
      	// execute application specific code
		task1_code();

		// it would be nice to check if we missed a deadline here... why don't
		// you try by yourself?

        struct timespec time_t1; // struct to store the current time in the characteristic function of the first thread: void *task1( void *), in order to check if some deadlines have been missed
	    clock_gettime(CLOCK_REALTIME, &time_t1); // store the current time of the day in the time_1 variable

        if((time_t1.tv_sec>next_arrival_time[0].tv_sec)||((time_t1.tv_sec==next_arrival_time[0].tv_sec)&&(time_t1.tv_nsec>next_arrival_time[0].tv_nsec))){
              missed_deadlines[0] = missed_deadlines[0]+1;
        }


		// sleep until the end of the current period (which is also the start of the
		// new one
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL); // make the thread sleep until a known given time (the beginning
		// of the next period: next_arrival_time[0])

		// the thread is ready and can compute the end of the current period for
		// the next iteration

        //just before exiting the loop, naturally I update the variable next_arrival_time[0] (the first loop is carried out, thanks to the fact that, before the loop, the beginning of the
        // second period is evaluated)
		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0]; // each time the thread is executed the subsequent arrival time (AKA beginning of the subsequent period)
		// is evaluated as seen before (N.B. before only the beginning of the second period was evaluated)
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task2_code()
{
	//print the id of the current task
  	printf(" 2[ "); fflush(stdout);
	int i,j;
	double uno;
  	for (i = 0; i < OUTERLOOP; i++)
    	{
      		for (j = 0; j < INNERLOOP; j++)
		{
			uno = rand()*rand()%10;
		}
    	}
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
		//same as task 1

		task2_code();

        struct timespec time_t2;
	    clock_gettime(CLOCK_REALTIME, &time_t2); // store the current time of the day in the time_1 variable

        if((time_t2.tv_sec>next_arrival_time[1].tv_sec)||((time_t2.tv_sec==next_arrival_time[1].tv_sec)&&(time_t2.tv_nsec>next_arrival_time[1].tv_nsec))){
              missed_deadlines[1] = missed_deadlines[1]+1;
        }


		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);

		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;

    	}
}

void task3_code()
{
	//print the id of the current task
  	printf(" 3[ "); fflush(stdout);
	int i,j;
	double uno;
  	for (i = 0; i < OUTERLOOP; i++)
    	{
      		for (j = 0; j < INNERLOOP; j++);		
			double uno = rand()*rand()%10;
    	}
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
		//same as task 1

		task3_code();

        struct timespec time_t3;
	    clock_gettime(CLOCK_REALTIME, &time_t3); // store the current time of the day in the time_1 variable

        if((time_t3.tv_sec>next_arrival_time[2].tv_sec)||((time_t3.tv_sec==next_arrival_time[2].tv_sec)&&(time_t3.tv_nsec>next_arrival_time[2].tv_nsec))){
              missed_deadlines[2] = missed_deadlines[2]+1;
        }

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);
		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
		next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
	
	}
}




//- APERIODIC TASKS

void task4_code()
{
  	printf(" 4[ "); fflush(stdout);
	for (int i = 0; i < OUTERLOOP; i++)
    	{
      		for (int j = 0; j < INNERLOOP; j++)
			double uno = rand()*rand();
    	}
  	printf(" ]4 "); fflush(stdout);
  	fflush(stdout);
}

void *task4( void *)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	//add an infinite loop 
	while (1)
    	{
		// wait for the proper condition to be signalled
		// See below why mutexes have been commented
//		pthread_mutex_lock(&mutex_task_4);
		pthread_cond_wait(&cond_task_4, &mutex_task_4);
//		pthread_mutex_unlock(&mutex_task_4);
		// execute the task code
 		task4_code();
	}
}

void task5_code()
{
  	printf(" 5[ "); fflush(stdout);
  	for (int i = 0; i < OUTERLOOP; i++)
    	{
      		for (int j = 0; j < INNERLOOP; j++)	
			double uno = rand()*rand();
    	}	
  	printf(" ]5 "); fflush(stdout);
}

void *task5( void *)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);
	//add an infinite loop  	
	while(1)
    	{
		// wait for the proper condition to be signalled
		// See below why mutexes have been commented
//      		pthread_mutex_lock(&mutex_task_5);
      		pthread_cond_wait(&cond_task_5, &mutex_task_5);
//      		pthread_mutex_unlock(&mutex_task_5);
		// execute the task code
      		task5_code();
    	}
}


void task6_code()
{
  	printf(" 6[ "); fflush(stdout);
	for (int i = 0; i < OUTERLOOP; i++)
    	{
      		for (int j = 0; j < INNERLOOP; j++)
			double uno = rand()*rand();
    	}
  	printf(" ]6 "); fflush(stdout);
  	fflush(stdout);
}

void *task6( void *)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	//add an infinite loop 
	while (1)
    	{
		// wait for the proper condition to be signalled
		// See below why mutexes have been commented
//		pthread_mutex_lock(&mutex_task_6);
		pthread_cond_wait(&cond_task_6, &mutex_task_6);
//		pthread_mutex_unlock(&mutex_task_6);
		// execute the task code
 		task6_code();
	}
}

