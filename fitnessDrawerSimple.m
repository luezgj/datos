exp_dir="/home/luezgj/Documentos/Datos/20-04-23/";

generations=1800;


fid = fopen ([exp_dir,"fitness.dat"],"r")
gen=1:generations;
best=zeros(1,generations);
prom=zeros(1,generations);

for generation=1:generations
	txt = fgetl (fid);
	txt=strrep(txt,"s","");
	vector=str2num(txt);
	best(generation)=vector(2);
	prom(generation)=vector(3);
endfor
                  
fclose(fid);
                                    
clf;
plot(gen,log2(prom),"r;prom;");
hold on;
plot(gen,log2(best),"b;best;");
yt = get(gca, 'YTick');
set (gca, 'YTickLabel', 2.^yt);	
print([exp_dir,"fitness.png"]);