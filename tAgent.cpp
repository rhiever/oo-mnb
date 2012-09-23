/*
 * tAgent.cpp
 *
 * This file is part of the Simon Memory Game project.
 *
 * Copyright 2012 Randal S. Olson.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include "tAgent.h"

tAgent::tAgent()
{
	nrPointingAtMe=1;
	ancestor = NULL;
	for(int i=0;i<maxNodes;i++)
    {
		states[i]=0;
		newStates[i]=0;
	}
	bestSteps=-1;
	ID=masterID;
	masterID++;
	saved=false;
	hmmus.clear();
	nrOfOffspring=0;
	retired=false;
	food=0;
    totalSteps=0;
    numStates = numInputs + numOutputs + 2;
#ifdef useANN
	ANN=new tANN;
#endif
}

tAgent::~tAgent()
{
	for (int i = 0; i < hmmus.size(); ++i)
    {
        delete hmmus[i];
    }
    
	if (ancestor!=NULL)
    {
		ancestor->nrPointingAtMe--;
		if (ancestor->nrPointingAtMe == 0)
        {
			delete ancestor;
        }
	}
#ifdef useANN
	delete ANN;
#endif
}

void tAgent::setupRandomAgent(int nucleotides)
{
	genome.resize(nucleotides);
    
	for(int i = 0; i < nucleotides; ++i)
    {
        bool createdHMGStartCodon = false;
        
        do
        {
            genome[i] = rand() & 255;
            
            createdHMGStartCodon = (i > 1 && genome[i - 1] == 42 && genome[i] == (255 - 42));

        } while (createdHMGStartCodon);
    }
    
	ampUpStartCodons();
    setupPhenotype();
#ifdef useANN
	ANN->setup();
#endif
}

void tAgent::loadAgent(char* filename)
{
	FILE *f=fopen(filename,"r");
	int i;
	genome.clear();
	while(!(feof(f)))
    {
		fscanf(f,"%i	",&i);
		genome.push_back((unsigned char)(i&255));
	}
    fclose(f);
	setupPhenotype();
}

void tAgent::loadAgentWithTrailer(char* filename)
{
#ifdef useANN
	ANN=new tANN;
	ANN->load(filename);
#else
	FILE *f=fopen(filename,"r+t");
	int i;
	genome.clear();
	fscanf(f,"%i	",&i);
	while(!(feof(f))){
		fscanf(f,"%i	",&i);
		genome.push_back((unsigned char)(i&255));
	}
	setupPhenotype();
#endif
}


void tAgent::ampUpStartCodons(void)
{
	int i,j;
    
#ifdef directedMutations
    // add start gates
    j = 0;
    
	for(i = 0; i < 4; ++i)
	{
		genome[j] = 42;
		genome[j + 1] = (255 - 42);
        //j+2 gate type
        //j+3 # ins
        //j+4 # outs
        for (int offset = 5; offset < 13; ++offset)
        {
            genome[j + offset] = rand() % numStates;
        }
        
        j += 270;
	}
    
    genome.resize(j);
#else
    // add start gates
	for(i = 0; i < 4; ++i)
	{
		j=rand()%((int)genome.size()-100);
		genome[j]=42;
		genome[j+1]=(255-42);
		for(int k=2;k<20;k++)
        {
			genome[j+k]=rand()&255;
        }
	}
    
#endif
}

void tAgent::inherit(tAgent *from, double mutationsPerInherit, double duplicationRate, double deletionRate, double addStateMutationRate, double removeStateMutationRate, int theTime)
{
	int nucleotides = (int)from->genome.size();
    int lastHMGStartCodonIndex = 0;
	double mutationRate = mutationsPerInherit / from->genome.size();
	vector<unsigned char> buffer;
	born = theTime;
    numHMGs = from->numHMGs;
    numStates = from->numStates;
	//ancestor=from;
	//from->nrPointingAtMe++;
	from->nrOfOffspring++;
	genome.clear();
	genome.resize(from->genome.size());
    
#ifdef directedMutations
    
    // gate per-site mutation
	for(int i = 0; i < nucleotides; ++i)
    {
        // disallow point mutations from deleting HMGs
        bool isHMGStartCodon = (from->genome[i] == 42 && from->genome[i + 1] == (255 - 42)) ||
                                (i > 1 && genome[i - 1] == 42 && from->genome[i] == (255 - 42));
        
        if (isHMGStartCodon)
        {
            lastHMGStartCodonIndex = i;
        }
        
		if(!isHMGStartCodon && randDouble < mutationRate)
        {
            bool createdHMGStartCodon = false;
            
            // disallow point mutations from creating new HMGs
            do
            {
                if ( (i - lastHMGStartCodonIndex) <= 11 && (i - lastHMGStartCodonIndex) >= 4 )
                {
                    genome[i] = rand() % numStates;
                }
                else
                {
                    genome[i] = rand() & 255;
                }
                
                createdHMGStartCodon = (genome[i] == 42 && from->genome[i + 1] == (255 - 42)) ||
                                        (i > 1 && genome[i - 1] == 42 && genome[i] == (255 - 42));

            } while (createdHMGStartCodon);
        }
		else
        {
			genome[i] = from->genome[i];
        }
    }
    
    // add gate
    if(randDouble < duplicationRate)
    {
        ++numHMGs;
        
        // duplicate gate
        if (rand() % 2 == 0)
        {
            int geneToDuplicate = 1 + (rand() % numHMGs);
            int geneCount = 0;
            int copyStartIndex = 0, copyEndIndex = 0;
            
            for (int i = 0; i < nucleotides; ++i)
            {
                bool isHMG = genome[i] == 42 && genome[i + 1] == (255 - 42);
                
                if (isHMG)
                {
                    ++geneCount;
                    
                    if (geneCount == geneToDuplicate)
                    {
                        copyStartIndex = i;
                        
                        if (isHMG)
                        {
                            copyEndIndex = i + 270;
                        }
                        else
                        {
                            copyEndIndex = copyStartIndex;
                        }
                        
                        break;
                    }
                }
            }
            
            buffer.clear();
            buffer.resize(0);
            buffer.insert(buffer.begin(), genome.begin() + copyStartIndex, genome.begin() + copyEndIndex);
        }
        
        // random gate
        else
        {
            buffer.clear();
            buffer.resize(270);
            
            buffer[0] = 42;
            buffer[1] = 255 - 42;
            
            for (int i = 2; i < buffer.size(); ++i)
            {
                if (i >=5 && i <= 12)
                {
                    buffer[i] = rand() % numStates;
                }
                else
                {
                    buffer[i] = rand() & 255;
                }
            }
        }
        
        genome.insert(genome.end(), buffer.begin(), buffer.end());
    }
    
    // delete gate
    if(numHMGs > 1 && randDouble < deletionRate)
    {
        --numHMGs;
        int geneToDelete = 1 + (rand() % numHMGs);
        int geneCount = 0;
        int deleteStartIndex = 0, deleteEndIndex = 0;
        
        for (int i = 0; i < nucleotides; ++i)
        {
            bool isHMG = genome[i] == 42 && genome[i + 1] == (255 - 42);
            
            if (isHMG)
            {
                ++geneCount;
                
                if (geneCount == geneToDelete)
                {
                    deleteStartIndex = i;
                    
                    if (isHMG)
                    {
                        deleteEndIndex = i + 270;
                    }
                    else
                    {
                        deleteEndIndex = deleteStartIndex;
                    }
                    
                    break;
                }
            }
        }
        
        genome.erase(genome.begin() + deleteStartIndex, genome.begin() + deleteEndIndex);
    }
    
    // add state
    if (randDouble < addStateMutationRate)
    {
        ++numStates;
    }
    
    // remove state
    if (numStates > (numInputs + numOutputs) && randDouble < removeStateMutationRate)
    {
        --numStates;
    }
    
#else
    
    // per-site mutation
	for(int i = 0; i < nucleotides; ++i)
    {
		if(randDouble < mutationRate)
        {
			genome[i] = rand() & 255;
        }
		else
        {
			genome[i] = from->genome[i];
        }
    }
    
    // duplication
    if((randDouble < duplicationRate) && (genome.size() < 20000))
    {
        int copyEndIndex = 15 + rand() & 511;
        int copyStartIndex = rand() % ((int)genome.size() - copyEndIndex);
        int insertionIndex = rand() % (int)genome.size();
        buffer.clear();
        buffer.insert(buffer.begin(), genome.begin() + copyStartIndex, genome.begin() + copyStartIndex + copyEndIndex);
        genome.insert(genome.begin() + insertionIndex, buffer.begin(), buffer.end());
    }
    
    // deletion
    if((randDouble < deletionRate) && (genome.size() > 1000))
    {
        int deleteEndIndex = 15 + rand() & 511;
        int deleteStartIndex = rand() % ((int)genome.size() - deleteEndIndex);
        genome.erase(genome.begin() + deleteStartIndex, genome.begin() + deleteStartIndex + deleteEndIndex);
    }
    
#endif

	setupPhenotype();
	fitness = 0.0;
#ifdef useANN
	ANN->inherit(ancestor->ANN, mutationRate);
#endif
}

void tAgent::setupPhenotype(void)
{
	tHMMU *hmmu;
    
    numHMGs = 0;
    
	if(hmmus.size() != 0)
    {
		for(int i = 0; i < hmmus.size(); ++i)
        {
			delete hmmus[i];
        }
    }
    
	hmmus.clear();
    
	for(int i = 0; i < genome.size() - 1; ++i)
    {
		if((genome[i] == 42) && (genome[i + 1] == (255 - 42)))
        {
			hmmu=new tHMMU;
            
            // deterministic gate
            if (genome[i + 2] % 2 == 0)
            {
                hmmu->setupQuick(genome, i, numStates);
            }
			
            // probabilistic gate
			else
            {
                hmmu->setup(genome, i, numStates);
            }
            
			hmmus.push_back(hmmu);
            
            ++numHMGs;
		}
	}
}


void tAgent::retire(void)
{
	retired=true;
}

unsigned char * tAgent::getStatesPointer(void)
{
	return states;
}

void tAgent::resetBrain(void)
{
	for(int i=0;i<maxNodes;i++)
    {
		states[i]=0;
    }
#ifdef useANN
	ANN->resetBrain();
#endif
}

void tAgent::updateStates(void)
{
	for(vector<tHMMU*>::iterator it = hmmus.begin(), end = hmmus.end(); it != end; ++it)
    {
		(*it)->update(&states[0],&newStates[0]);
    }
    
	for(int i=0;i<maxNodes;i++)
    {
		states[i]=newStates[i];
		newStates[i]=0;
	}
	++totalSteps;
}

void tAgent::showBrain(void)
{
	for(int i=0;i<maxNodes;i++)
    {
		cout<<(int)states[i];
    }
	cout<<endl;
}

void tAgent::initialize(int x, int y, int d)
{
	//int i,j;
	//unsigned char dummy;
	xPos=x;
	yPos=y;
	direction=d;
	steps=0;
	/*
	if((rand()&1)==0){
		scramble[1]=2;
		scramble[2]=1;
	}
	*/
}

tAgent* tAgent::findLMRCA(void)
{
	tAgent *r,*d;
	if(ancestor==NULL)
		return NULL;
	else{
		r=ancestor;
		d=NULL;
		while(r->ancestor!=NULL){
			if(r->ancestor->nrPointingAtMe!=1)
				d=r;
			r=r->ancestor;
		}
		return d;
	}
}

void tAgent::saveFromLMRCAtoNULL(FILE *statsFile,FILE *genomeFile){
	if(ancestor!=NULL)
		ancestor->saveFromLMRCAtoNULL(statsFile,genomeFile);
	if(!saved){ 
		fprintf(statsFile,"%i	%i	%i	%f	%i	%f	%i	%i\n",ID,born,(int)genome.size(),fitness,bestSteps,(float)totalSteps/(float)nrOfOffspring,correct,incorrect);
		fprintf(genomeFile,"%i	",ID);
		for(int i=0;i<genome.size();i++)
			fprintf(genomeFile,"	%i",genome[i]);
		fprintf(genomeFile,"\n");
		saved=true;
	}
	if((saved)&&(retired)) genome.clear();
}

/*
void tAgent::saveLOD(FILE *statsFile,FILE *genomeFile){
	if(ancestor!=NULL)
		ancestor->saveLOD(statsFile,genomeFile);
#ifdef useANN
	fprintf(genomeFile,"%i	",ID);
	fprintf(statsFile,"%i	%i	%i	%f	%i	%f	%i	%i\n",ID,born,(int)genome.size(),fitness,bestSteps,(float)totalSteps/(float)nrOfOffspring,correct,incorrect);
	ANN->saveLOD(genomeFile);
#else	
	fprintf(statsFile,"%i	%i	%i	%f	%i	%f	%i	%i\n",ID,born,(int)genome.size(),fitness,bestSteps,(float)totalSteps/(float)nrOfOffspring,correct,incorrect);
	fprintf(genomeFile,"%i	",ID);
	for(int i=0;i<genome.size();i++)
		fprintf(genomeFile,"	%i",genome[i]);
	fprintf(genomeFile,"\n");
#endif
	
}*/

void tAgent::showPhenotype(void)
{
	for(int i = 0; i < hmmus.size(); ++i)
    {
		hmmus[i]->show();
    }
	cout<<"------"<<endl;
}

void tAgent::saveToDot(const char *filename)
{
	FILE *f=fopen(filename,"w+t");
	int i,j,k,node;
	fprintf(f,"digraph brain {\n");
	fprintf(f,"	ranksep=2.0;\n");
    
    // determine which nodes to print (no connections = do not print)
    bool print_node[maxNodes];
    
    for(i = 0; i < maxNodes; i++)
    {
        print_node[i] = false;
    }
    
    for(i=0;i<hmmus.size();i++)
    {
        for(j=0;j<hmmus[i]->ins.size();j++)
        {
            print_node[hmmus[i]->ins[j]] = true;
        }
        
        for(k=0;k<hmmus[i]->outs.size();k++)
        {
            print_node[hmmus[i]->outs[k]] = true;
        }
    }
    
    for (int i = 0; i < numColors; ++i)
    {
        print_node[i] = true;
    }
    
    for (int i = maxNodes - 1; i >= (maxNodes - numOutputs); --i)
    {
        print_node[i] = true;
    }
    
    // input layer
	for(node=0;node<numColors;node++)
    {
        if(print_node[node])
        {
            fprintf(f,"	%i [shape=invtriangle,style=filled,color=cornflowerblue];\n",node);
        }
    }
    
    // hidden states
    for(node=numColors;node<maxNodes-numOutputs;node++)
    {
        if(print_node[node])
        {
            fprintf(f,"	%i [shape=circle,color=black];\n",node);
        }
    }
    
    // outputs
	for(node=maxNodes-numOutputs;node<maxNodes;node++)
    {
		fprintf(f,"	%i [shape=circle,style=filled,color=green];\n",node);
    }
    
    // connections
	for(i=0;i<hmmus.size();i++)
    {
		for(j=0;j<hmmus[i]->ins.size();j++)
        {
			for(k=0;k<hmmus[i]->outs.size();k++)
            {
				fprintf(f,"	%i	->	%i;\n",hmmus[i]->ins[j],hmmus[i]->outs[k]);
            }
		}
	}
    
    // which nodes go on the same level
    // inputs
    fprintf(f,"	{ rank=same; ");
    
    for(node = 0; node < numColors; node++)
    {
        if(print_node[node])
        {
            fprintf(f, "%d; ", node);
        }
    }
    
    fprintf(f, "}\n");
    
    // hidden states
    fprintf(f,"	{ rank=same; ");
    
    for(node = numColors; node < maxNodes-numOutputs; node++)
    {
        if(print_node[node])
        {
            fprintf(f, "%d; ", node);
        }
    }
        
    fprintf(f, "}\n");
    
    // outputs
    fprintf(f,"	{ rank=same; ");
    
    for(node = maxNodes-numOutputs; node < maxNodes; node++)
    {
        if(print_node[node])
        {
            fprintf(f, "%d; ", node);
        }
    }
    
    fprintf(f, "}\n");
	fclose(f);
}

void tAgent::saveToDotFullLayout(char *filename){
	FILE *f=fopen(filename,"w+t");
	int i,j,k;
	fprintf(f,"digraph brain {\n");
	fprintf(f,"	ranksep=2.0;\n");
	for(i=0;i<hmmus.size();i++){
		fprintf(f,"MM_%i [shape=box]\n",i);
		for(j=0;j<hmmus[i]->ins.size();j++)
			fprintf(f,"	t0_%i -> MM_%i\n",hmmus[i]->ins[j],i);
		for(k=0;k<hmmus[i]->outs.size();k++)
			fprintf(f,"	MM_%i -> t1_%i\n",i,hmmus[i]->outs[k]);
		
	}
	fprintf(f,"}\n");
}

void tAgent::setupDots(int x, int y,double spacing){
	double xo,yo;
	int i,j,k;
	xo=(double)(x-1)*spacing;
	xo=-(xo/2.0);
	yo=(double)(y-1)*spacing;
	yo=-(yo/2.0);
	dots.resize(x*y);
	k=0;
	for(i=0;i<x;i++)
		for(j=0;j<y;j++){
//			dots[k].xPos=(double)(rand()%(int)(spacing*x))+xo;
//			dots[k].yPos=(double)(rand()%(int)(spacing*y))+yo;
			dots[k].xPos=xo+((double)i*spacing);
			dots[k].yPos=yo+((double)j*spacing);
//			cout<<dots[k].xPos<<" "<<dots[k].yPos<<endl;
			k++;
		}
}

void tAgent::saveLogicTable(const char *filename)
{
    FILE *f=fopen(filename, "w");
	int i,j;
    
    fprintf(f,"s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,p15,,o1,o2\n");
    //fprintf(f,"s11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,,o1,o2\n");
    
    for(i = 0; i < (int)pow(2.0, 13.0); i++)
    {
        map<vector<int>, int> outputCounts;
        const int NUM_REPEATS = 1001;
        
        for (int repeat = 1; repeat < NUM_REPEATS; ++repeat)
        {
            for(j = 0; j < 30; j++)
            {
                if (j < 12)
                {
                    if(repeat == 1)
                    {
                        fprintf(f,"%i,",(i >> j) & 1);
                    }
                    
                    states[j] = (i >> j) & 1;
                }
                else if (j == 15)
                {
                    if(repeat == 1)
                    {
                        fprintf(f,"%i,",(i >> 12) & 1);
                    }
                    
                    states[j] = (i >> 12) & 1;
                }
                else
                {
                    states[j] = 0;
                }
            }
            
            updateStates();
            
            vector<int> output;
            // order: 30 31
            output.push_back(states[30]);
            output.push_back(states[31]);
            
            if (outputCounts.count(output) > 0)
            {
                outputCounts[output]++;
            }
            else
            {
                outputCounts[output] = 1;
            }
            
            // all repeats completed; determine most common output
            if (repeat == (NUM_REPEATS - 1))
            {
                map<vector<int>, int>::iterator it;
                map<vector<int>, int>::iterator mostCommonOutput = outputCounts.begin();
                
                for (it = outputCounts.begin(); it != outputCounts.end(); ++it)
                {
                    if (it->second > mostCommonOutput->second)
                    {
                        mostCommonOutput = it;
                    }
                }
                
                fprintf(f, ",%i,%i\n", mostCommonOutput->first[0], mostCommonOutput->first[1]);
            }
        }
	}
    
    fclose(f);
}

// saves the Markov network brain genome to a text file
void tAgent::saveGenome(const char *filename)
{
    FILE *f = fopen(filename, "w");
    
	for (int i = 0, end = (int)genome.size(); i < end; ++i)
    {
		fprintf(f, "%i	", genome[i]);
    }
    
	fprintf(f, "\n");
    
    fclose(f);
}
