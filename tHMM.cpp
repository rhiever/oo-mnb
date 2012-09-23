/*
 * tHMM.cpp
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

#include <stdlib.h>
#include "tHMM.h"

tHMMU::tHMMU(){
}

tHMMU::~tHMMU(){
	hmm.clear();
	sums.clear();
	ins.clear();
	outs.clear();
}
void tHMMU::setup(vector<unsigned char> &genome, int start, int numStates)
{
	int i,j,k;
	ins.clear();
	outs.clear();
	k=(start+3)%(int)genome.size();

	_xDim=1+(genome[(k++)%genome.size()]&3);
	_yDim=1+(genome[(k++)%genome.size()]&3);
	//cout<<"setup "<<(int)genome[start+2]<<" "<<(int)xDim<<" "<<(int)yDim<<endl;
	ins.resize(_yDim);
	outs.resize(_xDim);

	for(i=0;i<_yDim;++i)
    {
		ins[i]=genome[(k+i)%genome.size()] % numStates;
    }
    
	for(i=0;i<_xDim;++i)
    {
		outs[i]=genome[(k+4+i)%genome.size()] % numStates;
    }
	
	k += 8;
	hmm.resize(1<<_yDim);
	sums.resize(1<<_yDim);
    
	for(i=0;i<(1<<_yDim);++i)
    {
		hmm[i].resize(1<<_xDim);
        
		for(j=0;j<(1<<_xDim);++j)
        {
//			hmm[i][j]=(genome[(k+j+((1<<yDim)*i))%genome.size()]&1)*255;
			hmm[i][j]=genome[(k+j+((1<<_xDim)*i))%genome.size()];
            
			if(hmm[i][j]==0)
            {
                hmm[i][j]=1;
            }
            
			sums[i]+=hmm[i][j];
		}
	}
}

void tHMMU::setupQuick(vector<unsigned char> &genome, int start, int numStates)
{
	int i,j,k;
	ins.clear();
	outs.clear();
	k=(start+3) % (int)genome.size();
	
	_xDim=1+(genome[(k++)%genome.size()]&3);
	_yDim=1+(genome[(k++)%genome.size()]&3);
	//cout<<"setup "<<(int)genome[start+2]<<" "<<(int)xDim<<" "<<(int)yDim<<endl;

	ins.resize(_yDim);
	outs.resize(_xDim);
    
	for(i=0;i<_yDim;++i)
    {
		ins[i]=genome[(k+i)%genome.size()] % numStates;
    }

	for(i=0;i<_xDim;++i)
    {
		outs[i]=genome[(k+4+i)%genome.size()] % numStates;
    }
	
	k += 8;
	hmm.resize(1<<_yDim);
	sums.resize(1<<_yDim);
    
	for(i=0;i<(1<<_yDim);++i)
    {
        int highestIndex = 0;
        
		hmm[i].resize(1<<_xDim);
        
		for(j=0;j<(1<<_xDim);++j)
        {
			hmm[i][j] = 0;
            
            if ((int)genome[(k+j+((1<<_xDim)*i))%genome.size()] >= (int)genome[(k+highestIndex+((1<<_xDim)*i))%genome.size()])
            {
                highestIndex = j;
            }
        }
        
		//hmm[i][genome[(k+j+((1<<_xDim)*i))%genome.size()]&((1<<_xDim)-1)]=255;
        hmm[i][highestIndex] = 255;
		sums[i] = 255;
	}
}

void tHMMU::update(unsigned char *states, unsigned char *newStates)
{
	int I=0;
	int i,j,r;
    
	for(vector<int>::iterator it = ins.begin(), end = ins.end(); it != end; ++it)
    {
		I=(I<<1)+((states[*it])&1);
    }
    
	r=1+(rand()%(sums[I]-1));
	j=0;
    //	cout<<I<<" "<<(int)hmm.size()<<" "<<(int)hmm[0].size()<<endl;
	while(r > hmm[I][j])
    {
		r -= hmm[I][j];
		++j;
	}
    
	for(i = 0; i < outs.size(); ++i)
    {
		newStates[outs[i]] |= (j >> i) & 1;
    }
}

void tHMMU::show()
{
	int i,j;
	cout<<"INS: ";
	for(i=0;i<ins.size();i++)
		cout<<(int)ins[i]<<" ";
	cout<<endl;
	cout<<"OUTS: ";
	for(i=0;i<outs.size();i++)
		cout<<(int)outs[i]<<" ";
	cout<<endl;
	for(i=0;i<hmm.size();i++){
		for(j=0;j<hmm[i].size();j++)
			cout<<" "<<(double)hmm[i][j]/sums[i];
		cout<<endl;
	}
	cout<<endl;

/*
	for(i=0;i<hmm.size();i++){
		for(j=0;j<hmm[i].size();j++)
			cout<<(int)hmm[i][j]<<" ";
		cout<<endl;
	}
	*/
//	cout<<"------"<<endl;
}
