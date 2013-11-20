%parameters are: a_min, a_max, d, r, tau_c, c_1
%N_tilda_A, nu_1_A, nu_2_A, nu_3_A, nu_4_A (functional assumptions)
%N_tilda_Omega, nu_1_Omega, nu_2_Omega, nu_3_Omenga, nu_4_Omega (part 2)
%Total of 16 parameters
%
%global dynamic variables are: N_A, N_Omega
%agent dynamic variables are: a, m, c, agent_types

num_agent_types = 2;

%for now, set them all the same

a_min = 0.002;
a_max = 1.0;
d = 1.05;
r = 1.1;
tau_c = 48;
c_1 = 8;
N_tilda_A = 10^5;
N_tilda_Omega = 10^5;

%more intuitive functional form variables
f_zero_A = 0.5;
f_half_A = 0.45;
f_full_A = 0.05;
f_inf_A = 0;

f_zero_Omega = 0.5;
f_half_Omega = 0.3;
f_full_Omega = 0.1;
f_inf_Omega = 0;


%convert intuitive types to nu
h_one_A = 1/(f_zero_A - f_half_A);
h_two_A = 1/(f_half_A - f_inf_A);
h_three_A = 1/(f_full_A - f_inf_A);

h_one_Omega = 1/(f_zero_Omega - f_half_Omega);
h_two_Omega = 1/(f_half_Omega - f_inf_Omega);
h_three_Omega = 1/(f_full_Omega - f_inf_Omega);

nu_1_A = (h_one_A*h_three_A - h_two_A.^2)/(h_one_A + h_three_A - h_two_A);
nu_2_A = h_one_A - nu_1_A;
nu_3_A = log((h_three_A - nu_1_A)/nu_2_A); %remember that matlab uses natural log
nu_4_A = f_inf_A;

nu_1_Omega = (h_one_Omega*h_three_Omega - h_two_Omega.^2)/(h_one_Omega + h_three_Omega - h_two_Omega);
nu_2_Omega = h_one_Omega - nu_1_Omega;
nu_3_Omega = log((h_three_Omega - nu_1_Omega)/nu_2_Omega); %remember that matlab uses natural log
nu_4_Omega = f_inf_Omega;

%now we need to specify 16*num_agent_types of parameters

param_mat = zeros(16,num_agent_types);

param_mat(1,:) = a_min;
param_mat(2,:) = a_max;
param_mat(3,:) = d;
param_mat(4,:) = r;
param_mat(5,:) = tau_c;
param_mat(6,:) = c_1;

param_mat(7,:) = N_tilda_A;
param_mat(8,:) = nu_1_A;
param_mat(9,:) = nu_2_A;
param_mat(10,:) = nu_3_A;
param_mat(11,:) = nu_4_A;

param_mat(12,:) = N_tilda_Omega;
param_mat(13,:) = nu_1_Omega;
param_mat(14,:) = nu_2_Omega;
param_mat(15,:) = nu_3_Omega;
param_mat(16,:) = nu_4_Omega;

trans_fun = @(p,x) p(5) + 1./(p(2) + p(3)*exp(p(4)*x./p(1)));

%some initial parameters
%ini_N_A = 10^4;
%ini_N_Omega = 10^4;
ini_N = 10^3;
max_t = 10^2;
max_N = 10^4;

%global dynamic variables
N_A = zeros(1,num_agent_types);
N_Omega = zeros(1,num_agent_types);

%empty local dynamic variables
free_indexes = zeros(max_N,1);
last_free_index = 0;

a = zeros(max_N,1);
m = zeros(max_N,1);
c = zeros(max_N,1);
agent_types = zeros(max_N,1);

%initialize the population

free_indexes(1:(max_N - 4*ini_N)) = (4*ini_N + 1):max_N;
last_free_index = (max_N - 4*ini_N);
a(1:4*ini_N) = a_min + (a_max - a_min)*rand(4*ini_N,1);
c(1:4*ini_N) = c_1;

agent_types(1:2*ini_N) = 1;
m(1:ini_N) = 1;
m((ini_N + 1):2*ini_N) = 2;
N_A(1) = ini_N;
N_Omega(1) = ini_N;

agent_types((2*ini_N+1):4*ini_N) = 1;
m((2*ini_N+1):3*ini_N) = 1;
m((3*ini_N+1):4*ini_N) = 2;
N_A(2) = ini_N;
N_Omega(2) = ini_N;

tic;
for t = 1:max_t,
	for agent_num = randperm(max_N), %randomize ordering
		a_type = agent_types(agent_num);
		if (m(agent_num) == 1),
			omega = (param_mat(1,a_type)./a(agent_num))*trans_fun(param_mat(7:11,a_type),N_Omega);
			if rand < omega, %change agent microenvironment
				m(agent_num) = 2;
				N_A(a_type) = N_A(a_type) - 1;
				N_Omega(a_type) = N_Omega(a_type) + 1;
				c(agent_num) = param_mat(6,a_type);
			else,
				a(agent_num) = min(param_mat(4,a_type)*a(agent_num),param_mat(2,a_type));
			end;
		elseif (m(agent_num) == 2),
			alpha = (a(agent_num)./param_mat(2,a_type))*trans_fun(param_mat(12:16,a_type),N_A);
			if ((c(agent_num) < param_mat(6,a_type)) & (rand < alpha)),
				%change agent microenvironment
				m(agent_num) = 1;
				N_A(a_type) = N_A(a_type) + 1;
				N_Omega(a_type) = N_Omega(a_type) - 1;
			else,
				if a(t) > a_min,
					c(agent_num) = c(agent_num) + 1;
					a(agent_num) = a(agent_num)./param_mat(3,a_type);
					if (c(agent_num) > param_mat(5,a_type)), %reproduce?
						c(agent_num) = 0;
						
						%begin cell duplication
						new_agent_num = free_indexes(last_free_index);
						last_free_index = last_free_index - 1;
						
						agent_types(new_agent_num) = a_type;
						m(new_agent_num) = 2;
						a(new_agent_num) = param_mat(1,a_type) + (param_mat(2,a_type) - param_mat(1,a_type))*rand;
						c(new_agent_num) = 0;
					end;
				else, %cell differentiates
					free_indexes(last_free_index + 1) = agent_num;
					m(agent_num) = 0;
					last_free_index = last_free_index + 1;
				end;
			end;
		end;
	end;
end;
toc;
