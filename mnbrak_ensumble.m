%%% -----------------------------------------------------------------
% 'mnbrak.m' - this function performs initial bracketing of the minimum.
% It has been brutally adopted from the corresponding subroutine in OPUS.

function [AX, FA, BX, FB, CX, FC] = mnbrak_ensumble(AX, BX, FA, FB, f1, f2, lambda, df, which, MAXITER)
%%% -----------------------------------------------------------------
% INPUT:    AX, BX   -->  initial guesses for x;
%           FA, FB   -->  the corresponding function evaluations of AX and BX using FUNC
%           f1, f2   -->  the reference constitutive relations of both state variables at each iteration of
%                         optimal reconstruction algorithm (ODE system will be solved
%                         based on this constitutive relation)!
%           df       -->  the descent (ascent) direction:  which is the L2
%                         or H1 gradient of the cost function with respect to the
%                         constitutive relation, this descent direction is
%                         in the direction of f1 or f2
%          which     -->  an integer (1 or 2) indicating which of the
%                         constitutive relations we are considering! if 1,
%                         f1 is updated and line search in the f1 space is
%                         performed.
%          MAXITER   -->  Maximum number of bracketing iterations.
% OUTPUT: AX, FA, BX, FB, CX, FC --> the bracketing triplet and their corresponding function evaluations
% ------------------------------------------------------------------------
%                        Bracketing Routine
% ------------------------------------------------------------------------
% The calling convention of this function is slightly modified from that
% in Numerical Recipies. Also, it assumes that the minimum is in the direction
% of BX from AX, so the routine will not explore to the other side of BX!
% Note that this gets caught in an endles loop if this direction is uphill!
% INPUT:  AX, BX are initial guesses of X which might bracket the minimum,
%   FA, FB must be provided as the initial function values at points AX, BX.
%/////////////////////////////////////////////////////////////////////////
%   Update made by Bartek Protas on the 02nd Feb. 2001 
%   MAXITER - maximum allowed number of iterations; If the maximum number
%   of iterations is exceeded without bracketing the minimum, the function 
%   returns CX=BX=AX;


GOLD=1.618034; 
GLIMIT=100.0; 
TINY=1.E-20;
CX=0.0;
GOLDINV = 1.0/(1.0+GOLD);
ICOUNT=0;

while (FB > FA) 
  if (ICOUNT <= 8)
    %disp('      Taking a smaller step in direction of gradient.')
    CX = BX;
    FC = FB;
    BX = GOLDINV * (CX + GOLD*AX);
    FB = Jeval_ensumble(BX, f1, f2, lambda, df, which); % evaluate cost function
    ICOUNT = ICOUNT + 1;
  else
    %disp('Oops - it appears as if I am proceeding uphill.')
    %disp('I will now use a safer bracketing procedure.')
    DUM = AX;
    AX = BX;
    BX = DUM;
    DUM = FB;
    FB = FA;
    FA = DUM;
    CX = 0.0;
  end
end

if (CX == 0.0)
%The following factor of GOLD was reduced to (GOLD-1.) to save one function
%evaluation near convergence.
  CX = BX + (GOLD-1.0) * (BX-AX);
  FC = Jeval_ensumble(CX, f1, f2, lambda, df, which);

%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
  ICOUNT = ICOUNT + 1;
  if ( ICOUNT >= MAXITER ) 
    BX = AX;
    CX = AX;
    %disp('Cound not bracket the minimum! Bugging out ...')
    return;
  end
%%/////////////////////////////////////////////////////////////////////////
end

% here label = 1;
while ( FB >= FC ) 
  goto = 0;
  R = (BX-AX)*(FB-FC);
  Q = (BX-CX)*(FB-FA);
  U = BX-((BX-CX)*Q-(BX-AX)*R) / (2.*max(abs(Q-R),TINY)*sign(Q-R));
  ULIM = BX+GLIMIT*(CX-BX);
% this is stupid, but I couldn't figure out a different way to implement
% jumps out an IF-block.
  for dd=1:1; if ( ((BX-U) * (U-CX)) > 0.0 )
    FU = Jeval_ensumble(U, f1, f2, lambda, df, which);
    
%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
    ICOUNT=ICOUNT+1;
    if ( ICOUNT >= MAXITER ) 
      BX = AX;
      CX = AX;
      %disp('Cound not bracket the minimum! Bugging out ...')
      return;
    end
%%/////////////////////////////////////////////////////////////////////////
    goto = 0;
    if ( FU < FC )
      AX = BX;
      FA = FB;
      BX = U;
      FB = FU;
      goto = 1;
    elseif ( FU > FB )
      CX = U;
      FC = FU;
      goto = 1;
    end
    if ( goto == 1 ) break; end

    U = CX+GOLD*(CX-BX);
    FU = Jeval_ensumble(U, f1, f2, lambda, df, which);


%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
    ICOUNT = ICOUNT+1;
    if ( ICOUNT >= MAXITER ) 
      BX = AX;
      CX = AX;
      %disp('Cound not bracket the minimum! Bugging out ...')
      return;
    end
%%/////////////////////////////////////////////////////////////////////////
  elseif ( (CX-U)*(U-ULIM) > 0.0 )
    FU = Jeval_ensumble(U, f1, f2, lambda, df, which);

%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
    ICOUNT = ICOUNT+1;
    if ( ICOUNT >= MAXITER )
      BX = AX;
      CX = AX;
      %disp('Cound not bracket the minimum! Bugging out ...')
      return;
    end
%%/////////////////////////////////////////////////////////////////////////
    if ( FU < FC )
      BX = CX;
      CX = U;
      U = CX+GOLD*(CX-BX);
      FB = FC;
      FC = FU;
      FU = Jeval_ensumble(U, f1, f2, lambda, df, which);

%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
      ICOUNT = ICOUNT+1;
      if ( ICOUNT >= MAXITER ) 
        BX = AX;
	CX = AX;
	%disp('Cound not bracket the minimum! Bugging out ..')
        return;
      end
%%/////////////////////////////////////////////////////////////////////////
    end
  elseif ( (U-ULIM)*(ULIM-CX) >= 0.0 )
    U = ULIM;
    FU = Jeval_ensumble(U, f1, f2, lambda, df, which);

%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
    ICOUNT = ICOUNT+1;
    if ( ICOUNT >= MAXITER )
      BX = AX;
      CX = AX;
      %disp('Cound not bracket the minimum! Bugging out ...')
      return;
    end
%%/////////////////////////////////////////////////////////////////////////
  else
    U = CX+GOLD*(CX-BX);
    FU = Jeval_ensumble(U, f1, f2, lambda, df, which);
%%/////////////////////////////////////////////////////////////////////////
%%///// Update made by Bartek Protas on the 02nd Feb. 2001 ////////////////
    ICOUNT = ICOUNT+1;
    if ( ICOUNT >= MAXITER ) 
      BX = AX;
      CX = AX;
      %disp('Cound not bracket the minimum! Bugging out ...')
      return;
    end
    %%/////////////////////////////////////////////////////////////////////////
  end; end
  if ( goto ~= 1 )  
    AX = BX;
    BX = CX;
    CX = U;
    FA = FB;
    FB = FC;
    FC = FU;
  end

end
return;
