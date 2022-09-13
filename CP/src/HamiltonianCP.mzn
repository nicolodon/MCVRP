include "gecode.mzn";
include "globals.mzn";

% Multi-Capacitated Vehicle Routing problem - MiniZinc

% Variables

int: courier;
int: items;
array[1..courier] of int: load;
array[1..items] of 0..max(load): weight;
array[1..items+1,1..items+1] of int: dist;

set of int: VEHICLE = 1..courier;
set of int: CUSTOMER = 1..items;
set of int: LOAD = 0..max(load);
set of int: NODES = 1..items+2*courier;
set of int: DEPOT_NODES = items+1..items+2*courier;
set of int: START_DEPOT_NODES = items+1..items+courier;
set of int: END_DEPOT_NODES = items+courier+1..items+2*courier;

array[NODES] of int: demand = [  % Adapting demand to giant tour representation
  if i <= items then 
    weight[i]
  else 
    0
  endif
| i in NODES]; 

array[NODES, NODES] of int: distance = array2d(NODES,NODES,[           
  if i<=items /\ j <= items then 
    dist[i,j]
  elseif i<=items /\ j>items then % depot-customer
    dist[items+1,i]
  elseif j<=items /\ i>items then % customer-depot
    dist[j,items+1]
  else 
    dist[items+1,items+1] % depot-depot
  endif
   | i,j in NODES ]);  

array[NODES] of var NODES: successor; 
array[NODES] of var NODES: predecessor; % for redundant constraints
array[NODES] of var VEHICLE: vehicle; % which vehicle visits which customer?
array[NODES] of var LOAD: intload; % load when arriving at node n in NODES
var int: tot; % objective: total distance traversed by the couriers

% Constraints

% Initialization constraint
% Predecessor of start nodes are end nodes
constraint redundant_constraint(
   forall(n in (items+2..items+courier)) (
     predecessor[n] = n + courier-1
   )
);

constraint redundant_constraint(
   predecessor[items+1] = items+2*courier
);

% Successors of end nodes are start nodes
constraint 
   forall(n in (items+courier+1..items+2*courier-1)) (
     successor[n] = n-courier+1 
   );
constraint
   successor[items+2*courier] = items+1;

% Associate each start/end nodes with a vehicle
constraint 
   forall(n in START_DEPOT_NODES) (
     vehicle[n] = n-items
   );
   
constraint 
   forall(n in END_DEPOT_NODES) (
     vehicle[n] = n-items-courier
   );


% Vehicle load when starting at the depot
constraint 
   forall(n in START_DEPOT_NODES) (
     intload[n] = 0
   );
   

% Predecessor/successor constraints

constraint redundant_constraint(
   forall(n in NODES) (
      successor[predecessor[n]] = n
   )
);

constraint redundant_constraint(
   forall(n in NODES) (
      predecessor[successor[n]] = n
   )
);

% Alldiff + subtour elimination constraints
constraint 
   circuit(successor);
   
constraint redundant_constraint(
   circuit(predecessor)
);


% Vehicle constraints 

% Vehicle of node i is the same as the vehicle for the predecessor
constraint redundant_constraint(
   forall(n in CUSTOMER) (
      vehicle[predecessor[n]] = vehicle[n]
   )
);
constraint 
   forall(n in CUSTOMER) (
      vehicle[successor[n]] = vehicle[n]
   );

% Load constraints

constraint 
   forall(n in CUSTOMER) (
      intload[n] + demand[n] = intload[successor[n]]
   );
constraint
   forall(n in START_DEPOT_NODES) (
      intload[n] = intload[successor[n]] 
   );

% Check that partial load is always less or equal then the capacity for each vehicle.
constraint redundant_constraint(
    forall(i in CUSTOMER)(
        intload[i] <= load[vehicle[i]]
    )
);

% Check that final load is less or equal then the capacity for each vehicle.
constraint forall(i in VEHICLE)(
    intload[i+items+courier] <= load[i]
);


% Symmetry breaking constraints

% function var int: equality(var int: a, int: a1) = if a=a1 then 1 else 0 endif;
% constraint symmetry_breaking_constraint(forall(i,j in VEHICLE)(
%      (load[i]==load[j] /\ i>j) -> lex_less([equality(vehicle[k],i)| k in CUSTOMER], [equality(vehicle[w],j)| w in CUSTOMER]) 
%    )
% );


% function var int: samecouriersucc(var int: i, int: c) = if vehicle[i]=c then successor[i] else 0 endif;
% function var int: samecourierpred(var int: i, int: c) = if vehicle[i]=c then predecessor[i] else 0 endif;
% constraint symmetry_breaking_constraint(forall(c in VEHICLE)(
%       lex_less( [samecourierpred(j,c) | j in NODES],[samecouriersucc(i,c) | i in NODES]))
%  );

% Objective

constraint
 tot = sum([distance[i, successor[i]] |i in NODES]);

solve 
:: seq_search([int_search(successor, dom_w_deg, indomain_random),
               int_search(predecessor, dom_w_deg, indomain_random),
               int_search(vehicle, first_fail, indomain_random),
               int_search(intload, first_fail, indomain_min),
               relax_and_reconstruct(successor, 92), % 89, 92 have been used 
               restart_luby(items*items)
              ])   
minimize tot; % distance


% Output

output 
   [ "objective = "] ++ [show(tot)] ++
   [ ";\nvehicle = " ] ++ [ show(vehicle) ]++
   [ ";\nsuccessor = "] ++          [ show(successor) ] ++
   [ ";\npredecessor = "] ++ [show(predecessor)] ++
   [ ";\nload = "] ++ [show(intload)] ++  
   [ ";\n"]
;