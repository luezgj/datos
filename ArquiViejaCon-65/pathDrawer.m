exp_dir="/home/luezgj/Documentos/ArquiViejaCon-65/"

fid = fopen ([exp_dir,"bestCateg.dat"],"r")
clf;
hold on;
allx=cell(1,48);
ally=cell(1,48);

#y=objectY
#x=objectX-agentX

for trial=1:24
	y=zeros(1,884);
	x=zeros(1,884);
	for t=1:884
		txt = fgetl (fid);
		vector=str2num(txt);
		y(t)=vector(4);
		x(t)=vector(3)-vector(5);
	endfor
	ally{trial}=y;
	allx{trial}=x;
	plot(x,y,"r");
endfor

for trial=1:24
	y=zeros(1,884);
	x=zeros(1,884);
	for t=1:884
		txt = fgetl (fid);
		vector=str2num(txt);
		y(t)=vector(4);
		x(t)=vector(3)-vector(5);
	endfor
	ally{trial+24}=y;
	allx{trial+24}=x;
	plot(x,y,"b");
endfor

print([exp_dir,"agentWay.jpg"]);
fclose (fid)

for trial=1:24
	clf;
	hold on;
	plot(allx{trial},ally{trial},"r");
	plot(allx{trial+24},ally{trial+24},"b");
	print([exp_dir,"agentWay_",num2str(trial),".jpg"]);
endfor