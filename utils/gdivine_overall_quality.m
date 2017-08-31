function qual = gdivine_overall_quality(r)

% Function to compute overall quality given feature vector 'r'

load svm_params_Live.mat

%% Classification
atrain = repmat(a_gdivine,[size(r,1) 1]);btrain = repmat(b_gdivine,[size(r,1) 1]);
x_curr = atrain.*r+btrain;
[~,~,p] = svmpredict(1, x_curr, model_class_gdivine, '-b 1'); % 384 X 5 matrix
p(:,model_class_gdivine.Label) = p;

%% Regression
q=[];
for i=1:5
    [q(:,i), ~,~] = svmpredict(1,x_curr,model_reg_gdivine{i});
end

%% Final GGSM-Divine score
qual = sum(p.*q,2);
