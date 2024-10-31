%%% -----------------------------------------------------------------
% 'Brent.m' - this function performs the line search minimization using the Brent's method;
% ref: algorithms for minimization without derivatives by Brent.
% It has been brutally adopted from the corresponding subroutine in OPUS.

function [BX, FB] = Brent_ensumble_C3_3C_FullModel(AX, BX, CX, FA, FB, FC, TOL, ITMAX, f1, f2, lambda, df, which)
%%% -----------------------------------------------------------------
% INPUT
% AX, FA, BX, FB, CX, FC - the bracketing triplet and their cost function values;
% TOL                    - tolerance;
% ITMAX                  - maximum number of iterations;
% f1, f2                 - the reference constitutive relations at each iteration of
%                           optimal reconstruction algorithm (ODE system will be solved
%                           based on this constitutive relation) for both
%                           state variables. f1 and f2 are constitutive
%                           relations as a function of each of the state
%                           variables.
% df                     - the descent (ascent) direction:  which is the L2
%                           or H1 gradient of the cost function with respect to the
%                           constitutive relation. The descent direction is
%                           with respect to one of the constitutive
%                           relations
% which                  - an integer (1 or 2) determining which
%                          constitutive relation we are considering.
% WARNING! - the function name is hardwired in the code;
% OUTPUT
% BX     - the minimizer;
% FB     - the corresponding value of the functional;
% min.f -- version:  2.03
%
% The following two functions achieve a line minimization of the 
% user supplied function FUNC, which is presumed to be a function of
% a singe variable (all appropriate parameters, i.c.'s, b.c.'s, etc. 
% for the computation of FUNC must be passed to it with common
% blocks or, if using Fortran 90, modules.)  The following set
% of routines are written in Fortran 77 and are taken directly from
% Numerical Recipes with only slight modification of the calling 
% convention, which saves a few calls to FUNC if used properly.
%
% To use these routines, first set up an initial bracketing guess:
%
% AX=X_GUESS_A    <- Set this as a parameter.
% BX=X_GUESS_B    <- Set this as a parameter.
% FA=FUNC(AX)
% FB=FUNC(BX)
% CALL MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
%
% Then, call the minimization routine:
%
% CALL BRENT (AX,BX,CX,FA,FB,FC,FUNC,TOL,ITMAX)
%
% INPUT:  AX, BX, and CX are a bracketing triplet of points from MNBRAK,
%   and FB must be provided as the initial function value at the point BX.
%   TOL*BX is the accuracy of the final answer.
%   ITMAX is maximum number of iterations.
% OUTPUT: FB is the minimum value of FUNC, at the minimum point BX.

CGOLD=0.3819660; 
ZEPS=1.0E-20; 
D=0.0;

% Storing the value corresponding to ALPHA=0.0; will need it below;
F0 = FA;
FW=min(FA, FC);
if (FW == FA)
  W = AX;
  V = CX;
  FV = FC;
else
  W = CX;
  V = AX;
  FV= FA;
end
X = BX;
FX = FB;
A = min( AX, CX);
B = min( AX, CX);

for ITER=1:ITMAX
  goto = 0;
  if ( ITER <= 2 )
     E = 2.0 * (B-A);
  end
  XM = 0.5 * (A+B);
  TOL1 = TOL * abs(X) + ZEPS;
  TOL2 = 2.0 * TOL1;
%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 03rd Apr. 2001 ////////////////
%        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
% Here we use the Wolfe's sufficient decrease condition combined with the original condition;
% The Wolfe's condition is temporarily turned off and the standard check is used;
%  if ( (FX <= (F0-CWolfe*X*GxP)) | (abs(X-XM) <= (TOL2-.5*(B-A))) ) 
  if ( abs(X-XM) <= (TOL2-.5*(B-A)) ) 
    goto = 3;
    break;
  end
%%/////////////////////////////////////////////////////////////////////////
% this is stupid, but I could not figure out a different way to implement
% jumps out an IF-block.
  for dd=1:1; if ( (abs(E) > TOL1) | (ITER <= 2) ) 
    R = (X-W)*(FX-FV);
    Q = (X-V)*(FX-FW);
    P = (X-V)*Q-(X-W)*R;
    Q = 2.*(Q-R);
    if ( Q > 0.0) 
      P = -P;
    end
    Q = abs(Q);
    ETEMP = E;
    E = D;
    if( (abs(P) >= abs(0.5*Q*ETEMP)) | (P <= Q*(A-X)) | (P >= Q*(B-X)) ) 
      goto = 1;
      break;
    end
    D = P/Q;
    U = X+D;
    if ( ((U-A) <= TOL2) | ((B-U) <= TOL2) ) 
      D = TOL1 * sign(XM-X);
    end
    goto = 2;
  end; end

% here is label = 1;
  if ( goto ~= 2 )
    if ( X >= XM )
      E = A-X;
    else
      E = B-X;
    end
    D = CGOLD * E;
  end

% here is label = 2;
  if ( abs(D) >= TOL1 )
    U = X+D;
  else
    U = X+TOL1 * sign(D);
  end   
  FU = Jeval_ensumble_C3_3C_FullModel(U, f1, f2, lambda, df, which);
  
  if (FU <= FX) 
    if (U >= X)
      A = X;
    else
      B = X;
    end
    V = W;
    FV = FW;
    W = X;
    FW = FX;
    X = U;
    FX = FU;
  else
    if (U < X)
      A = U;
    else
      B = U;
    end
    if ( (FU <= FW) | (W == X) ) 
      V = W;
      FV = FW;
      W = U;
      FW = FU;
    elseif ( (FU <= FV) | (V == X) | (V == W) ) 
        V = U;
	FV = FU;
    end
  end
end
if ( goto == 0 )
  text = sprintf('No convergence in %d iterations',ITMAX);
  disp(text);
end

BX = X;
FB = FX;
return;


