function P = createP(A,B,C,omega, alpha,opts)

    P = opts.At'*(2*alpha*vec(omega) - 2*opts.c);
    P = P - 2*alpha*B;
    P = P';
end

% function P = createP(A,B,C,omega, alpha,opts)
%     P = [];
%     for i=1:size(A,3)
%         P = [P trace(A(:,:,i)' * (2*alpha*omega - 2*C))];
%     end
%     P = P - 2*alpha*B';
% end