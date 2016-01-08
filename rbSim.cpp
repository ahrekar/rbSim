#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "ParticleGenerator.h"

int main(int argc, char* argv[]){
	//int n = stoi(argv[1]); 
	int n = 1;

	std::vector<Particle> particles;
	ParticleGenerator pg();

	for(int i =0; i < n; i++){
		particles.insert(particles.begin(),pg.generateParticle());
	}

	return 0;
}
