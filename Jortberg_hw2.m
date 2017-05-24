%Problem 3
theta=1;
x = [0:0.1:10];

prob = theta*exp(-theta*x);
figure; plot(x, prob); xlabel('x'); ylabel('p(x|theta)'); title('theta fixed at 1')

x=2;
theta=[0:0.01:5];
prob = zeros(length(theta),1);
for i =1:length(theta)
    prob(i) = theta(i)*exp(-theta(i)*x);
end
figure; plot(theta, prob); xlabel('theta'); ylabel('p(x|theta)'); title('x=2, 0<=theta<=5')

%Problem 6
theta = [0:0.01:1];
p1_s0 = 2*(1-theta);
p1_s1 = 2*theta;
figure; plot(theta, p1_s0, 'b'); hold on, plot(theta, p1_s1, 'r');
title('theta v p(theta|D)'); xlabel('theta'); ylabel('p(theta|D)'); legend('s=0', 's=1')
