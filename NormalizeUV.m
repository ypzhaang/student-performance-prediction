
function [U, V] = NormalizeUV(U, V, NormV, Norm)
U = U;
V = V;
%     K = size(U,2);
%     if Norm == 2
%         if NormV
%             norms = max(1e-15,sqrt(sum(V.^2,1)))';
%             V = V*spdiags(norms.^-1,0,K,K);
%             U = U*spdiags(norms,0,K,K);
%         else
%             norms = max(1e-15,sqrt(sum(U.^2,1)))';
%             U = U*spdiags(norms.^-1,0,K,K);
%             V = V*spdiags(norms,0,K,K);
%         end
%     else
%         if NormV
%             norms = max(1e-15,sum(abs(V),1))';
%             V = V*spdiags(norms.^-1,0,K,K);
%             U = U*spdiags(norms,0,K,K);
%         else
%             norms = max(1e-15,sum(abs(U),1))';
%             U = U*spdiags(norms.^-1,0,K,K);
%             V = V*spdiags(norms,0,K,K);
%         end
%     end

        