exp_dir="C:/Users/lucho/Documents/Tema1/Datos/myMix/GAParameterTest/";

generations=2000;

Tou=cell(1,5);
CrosM=cell(1,2);
CrosPoi=cell(1,4);
Elite=cell(1,4);
Pop=cell(1,3);
Mut=cell(1,3);
CrosPro=cell(1,4);
Survivors=cell(1,4);


function processFile(prefix)
	
endfunction

for tourn=0:4
	for crossover=0:1	
		for eliteFract=0:3
			for popSize=0:2
				for mutVar=0:2
					for crosProb=0:3
						for survPerc=0:3
                			if (crossover == 0) 
								for crossPoints=0:3
									prefix= ["T",num2str(tourn),"_CM",num2str(crossover),"_CP",num2str(crossPoints),"_EF",num2str(eliteFract),"_PS",num2str(popSize),"_MV",num2str(mutVar),"_CP",num2str(crosProb),"_SP",num2str(survPerc)];
									fid = fopen ([exp_dir,prefix,"fitness.dat"],"r")
									gen=zeros(1,generations);
									best=zeros(1,generations);
									prom=zeros(1,generations);

									for generation=1:generations
										txt = fgetl (fid);
										txt=strrep(txt,"s","");
										vector=str2num(txt);
										prom(generation)=vector(3);
										gen(generation)=generation;
										fgetl (fid);
										txt = fgetl (fid);
										vector=str2num(txt);
										best(generation)=vector(1);
									endfor
                  
                  fclose(fid);
                  
									endProm=mean(prom(end-9:end));

									Tou{tourn+1}=[Tou{tourn+1}, endProm];
									CrosM{crossover+1}=[CrosM{crossover+1}, endProm];
									CrosPoi{crossPoints+1}=[CrosPoi{crossPoints+1}, endProm];
									Elite{eliteFract+1}=[Elite{eliteFract+1}, endProm];
									Pop{popSize+1}=[Pop{popSize+1}, endProm];
									Mut{mutVar+1}=[Mut{mutVar+1}, endProm];
									CrosPro{crosProb+1}=[CrosPro{crosProb+1}, endProm];
									Survivors{survPerc+1}=[Survivors{survPerc+1}, endProm];
                  
                  
									clf;
									plot(gen,log2(prom),"r;prom;");
									hold on;
									plot(gen,log2(best),"b;best;");
									yt = get(gca, 'YTick');
									set (gca, 'YTickLabel', 2.^yt);	
									print([exp_dir,prefix,"fitness.jpg"]);
                  
                  
								endfor
							else
								prefix= ["T",num2str(tourn),"_CM",num2str(crossover),"_CP4_EF",num2str(eliteFract),"_PS",num2str(popSize),"_MV",num2str(mutVar),"_CP",num2str(crosProb),"_SP",num2str(survPerc)];
								fid = fopen ([exp_dir,prefix,"fitness.dat"],"r")
								gen=zeros(1,generations);
								best=zeros(1,generations);
								prom=zeros(1,generations);

								for generation=1:generations
									txt = fgetl (fid);
									txt=strrep(txt,"s","");
									vector=str2num(txt);
									prom(generation)=vector(3);
									gen(generation)=generation;
									fgetl (fid);
									txt = fgetl (fid);
									vector=str2num(txt);
									best(generation)=vector(1);
								endfor

                fclose(fid);
								endProm=mean(prom(end-9:end));

								Tou{tourn+1}=[Tou{tourn+1}, endProm];
								CrosM{crossover+1}=[CrosM{crossover+1}, endProm];
								CrosPoi{crossPoints+1}=[CrosPoi{crossPoints+1}, endProm];
								Elite{eliteFract+1}=[Elite{eliteFract+1}, endProm];
								Pop{popSize+1}=[Pop{popSize+1}, endProm];
								Mut{mutVar+1}=[Mut{mutVar+1}, endProm];
								CrosPro{crosProb+1}=[CrosPro{crosProb+1}, endProm];
								Survivors{survPerc+1}=[Survivors{survPerc+1}, endProm];

       
								clf;
								plot(gen,log2(prom),"r;prom;");
								hold on;
								plot(gen,log2(best),"b;best;");
								yt = get(gca, 'YTick');
								set (gca, 'YTickLabel', 2.^yt);	
								print([exp_dir,prefix,"fitness.jpg"]);
                
							endif
						endfor
					endfor
				endfor
			endfor
		endfor	
	endfor
endfor

pkg load communications;

tam=size(Tou{1})(2);
Tou=cell2mat(Tou);
Tou=vec2mat(Tou,tam);
Tou=transpose(Tou);

tam=size(CrosM{1})(2);
CrosM=cell2mat(CrosM);
CrosM=vec2mat(CrosM,tam);
CrosM=transpose(CrosM);

tam=size(CrosPoi{1})(2);
CrosPoi=cell2mat(CrosPoi);
CrosPoi=vec2mat(CrosPoi,tam);
CrosPoi=transpose(CrosPoi);

tam=size(Elite{1})(2);
Elite=cell2mat(Elite);
Elite=vec2mat(Elite,tam);
Elite=transpose(Elite);

tam=size(Pop{1})(2);
Pop=cell2mat(Pop);
Pop=vec2mat(Pop,tam);
Pop=transpose(Pop);

tam=size(Mut{1})(2);
Mut=cell2mat(Mut);
Mut=vec2mat(Mut,tam);
Mut=transpose(Mut);

tam=size(CrosPro{1})(2);
CrosPro=cell2mat(CrosPro);
CrosPro=vec2mat(CrosPro,tam);
CrosPro=transpose(CrosPro);

tam=size(Survivors{1})(2);
Survivors=cell2mat(Survivors);
Survivors=vec2mat(Survivors,tam);
Survivors=transpose(Survivors);

Tou=statistics (Tou);
CrosM=statistics(CrosM);
CrosPoi=statistics(CrosPoi);
Elite=statistics(Elite);
Pop=statistics(Pop);
Mut=statistics(Mut);
CrosPro=statistics(CrosPro);
Survivors=statistics(Survivors);

fid = fopen ([exp_dir,"staticsResume.txt"],"w")

fdisp (fid, "Tournament size:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, Tou);
fdisp (fid, "");

fdisp (fid, "Crossover mode:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, CrosM);
fdisp (fid, "");

fdisp (fid, "Crossover points:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, CrosPoi);
fdisp (fid, "");

fdisp (fid, "Elite fractions:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, Elite);
fdisp (fid, "");

fdisp (fid, "Population size:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, Pop);
fdisp (fid, "");

fdisp (fid, "Mutation variance:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, Mut);
fdisp (fid, "");

fdisp (fid, "Crossover probability:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, CrosPro);
fdisp (fid, "");

fdisp (fid, "Offspring survivors probability:");
fdisp (fid, "minimum \t fst quart \t median \t thrd quart \t maximum \t mean \t std dev \t skewness \t kurtosis");
fdisp (fid, Survivors);
fdisp (fid, "");

fclose (fid);