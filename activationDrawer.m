exp_dir="C:/Users/lucho/Documents/Tema1/Datos/myMix/primerPrueba/"

fid = fopen ([exp_dir,"bestCateg.dat"],"r")
clf;
hold off;
time=(0.1:0.1:88.4);
n1=cell(1,48);
n2=cell(1,48);
n3=cell(1,48);
n4=cell(1,48);
n5=cell(1,48);

#y=objectY
#x=objectX-agentX

for trial=1:24
	v1=zeros(1,884);
	v2=zeros(1,884);
	v3=zeros(1,884);
	v4=zeros(1,884);
	v5=zeros(1,884);
	for t=1:884
		txt = fgetl (fid);
		vector=str2num(txt);
		v1(t)=vector(13);
		v2(t)=vector(14);
		v3(t)=vector(15);
		v4(t)=vector(16);
		v5(t)=vector(17);
		endfor
	n1{trial}=v1;
	n2{trial}=v2;
	n3{trial}=v3;
	n4{trial}=v4;
	n5{trial}=v5;
	plot(time,v1,"r");
	print([exp_dir,"activation_",num2str(trial),"_c1.jpg"]);
	plot(time,v2,"r");
	print([exp_dir,"activation_",num2str(trial),"_c2.jpg"]);
	plot(time,v3,"r");
	print([exp_dir,"activation_",num2str(trial),"_c3.jpg"]);
	plot(time,v4,"r");
	print([exp_dir,"activation_",num2str(trial),"_c4.jpg"]);
	plot(time,v5,"r");
	print([exp_dir,"activation_",num2str(trial),"_c5.jpg"]);
endfor


for trial=1:24
	v1=zeros(1,884);
	v2=zeros(1,884);
	v3=zeros(1,884);
	v4=zeros(1,884);
	v5=zeros(1,884);
	for t=1:884
		txt = fgetl (fid);
		vector=str2num(txt);
		v1(t)=vector(13);
		v2(t)=vector(14);
		v3(t)=vector(15);
		v4(t)=vector(16);
		v5(t)=vector(17);
		endfor
	n1{trial+24}=v1;
	n2{trial+24}=v2;
	n3{trial+24}=v3;
	n4{trial+24}=v4;
	n5{trial+24}=v5;
	plot(time,v1,"b");
	print([exp_dir,"activation_",num2str(trial),"_l1.jpg"]);
	plot(time,v2,"b");
	print([exp_dir,"activation_",num2str(trial),"_l2.jpg"]);
	plot(time,v3,"b");
	print([exp_dir,"activation_",num2str(trial),"_l3.jpg"]);
	plot(time,v4,"b");
	print([exp_dir,"activation_",num2str(trial),"_l4.jpg"]);
	plot(time,v5,"b");
	print([exp_dir,"activation_",num2str(trial),"_l5.jpg"]);
endfor

fclose (fid)

for trial=1:24
	clf;
	hold on;
	plot(time,n1{trial},"r");
	plot(time,n1{trial+24},"b");
	print([exp_dir,"activation_",num2str(trial),"_1.jpg"]);
	clf;
	hold on;
	plot(time,n2{trial},"r");
	plot(time,n2{trial+24},"b");
	print([exp_dir,"activation_",num2str(trial),"_2.jpg"]);
	clf;
	hold on;
	plot(time,n3{trial},"r");
	plot(time,n3{trial+24},"b");
	print([exp_dir,"activation_",num2str(trial),"_3.jpg"]);
	clf;
	hold on;
	plot(time,n4{trial},"r");
	plot(time,n4{trial+24},"b");
	print([exp_dir,"activation_",num2str(trial),"_4.jpg"]);
	clf;
	hold on;
	plot(time,n5{trial},"r");
	plot(time,n5{trial+24},"b");
	print([exp_dir,"activation_",num2str(trial),"_5.jpg"]);
endfor