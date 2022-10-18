function data = legendreTable(points)
%% LegendreTable
%
% Gets the roots of the legendre polynomials and the weights for gaussian
% integration. The function input is the number of points and the output is
% a matriz with each row being for each point i and the first colum beeng
% the root and the second the weight. The legendre roots are defined in the
% interval [-1, 1], beeing always symetric about the origin. The function
% returns the roots and weights until a maximum number of points of 10.

switch points
    case 1
        data = [
        0.000000000000000 1.000000000000000];
    case 2
        data = [
         0.577350269189626 1.000000000000000
        -0.577350269189626 1.000000000000000];
    case 3
        data = [
         0.000000000000000 0.888888888888889
         0.774596669241483 0.555555555555556
        -0.774596669241483 0.555555555555556];
    case 4
        data = [
         0.339981043584856 0.652145154862546
         0.861136311594053 0.347854845137454
        -0.339981043584856 0.652145154862546
        -0.861136311594053 0.347854845137454];
    case 5
        data = [
         0.000000000000000 0.568888888888889 
         0.538469310105683 0.478628670499366 
         0.90617984593866  0.236926885056189
        -0.538469310105683 0.478628670499366 
        -0.90617984593866  0.236926885056189];
    case 6
        data = [
         0.238619186083197 0.467913934572691 
         0.661209386466265 0.360761573048139 
         0.932469514203152 0.171324492379170
        -0.238619186083197 0.467913934572691 
        -0.661209386466265 0.360761573048139 
        -0.932469514203152 0.171324492379170];
    case 7
        data = [
         0.000000000000000 0.417959183673469 
         0.405845151377397 0.381830050505119 
         0.741531185599394 0.279705391489277 
         0.949107912342759 0.129484966168870
        -0.405845151377397 0.381830050505119 
        -0.741531185599394 0.279705391489277 
        -0.949107912342759 0.129484966168870];
    case 8
        data = [
         0.183434642495650 0.362683783378362 
         0.525532409916329 0.313706645877887 
         0.796666477413627 0.222381034453374 
         0.960289856497536 0.101228536290376
        -0.183434642495650 0.362683783378362 
        -0.525532409916329 0.313706645877887 
        -0.796666477413627 0.222381034453374 
        -0.960289856497536 0.101228536290376];
    case 9
        data = [
         0.000000000000000 0.330239355001260 
         0.324253423403809 0.312347077040003 
         0.613371432700590 0.260610696402935 
         0.836031107326636 0.180648160694857 
         0.968160239507626 0.081274388361574
        -0.324253423403809 0.312347077040003 
        -0.613371432700590 0.260610696402935 
        -0.836031107326636 0.180648160694857 
        -0.968160239507626 0.081274388361574];
    case 10
        data = [
         0.148874338981631 0.295524224714753 
         0.433395394129247 0.269266719309996 
         0.679409568299024 0.219086362515982 
         0.865063366688985 0.149451349150581 
         0.973906528517172 0.066671344308688
        -0.148874338981631 0.295524224714753 
        -0.433395394129247 0.269266719309996 
        -0.679409568299024 0.219086362515982 
        -0.865063366688985 0.149451349150581 
        -0.973906528517172 0.066671344308688];
end