function output = loadNuisance(subNum, runNum)
% Find and import the nuisance-regressor file for this run
pths = specifyPaths;
subID = sprintf('sub-%02.f', subNum);
fname = findSubData(pths, subID, runNum, 'desc-confounds_timeseries.tsv');
nuis = readtable(fname, 'FileType', 'Text', 'Delimiter', '\t');
% There are a variable number of columns (e.g. multiple derivatives)
% And I really just care about a handful, so subset here.
% There are others like csf_wm, rot_z etc that MAY be useful? but idk
output = [nuis.global_signal, nuis.csf, nuis.white_matter];
% Also include the motion parameters and their derivatives:
% But scale each down to 0:1 based on data type
maxRot = max([nuis.rot_x, nuis.rot_y, nuis.rot_z], [], 'all');
maxRotD1 = max([nuis.rot_x_derivative1, nuis.rot_y_derivative1, nuis.rot_z_derivative1], [], 'all');
maxTrans = max([nuis.trans_x, nuis.trans_y, nuis.trans_z], [], 'all');
maxTransD1 = max([nuis.trans_x_derivative1, nuis.trans_y_derivative1, nuis.trans_z_derivative1], [], 'all');
output = [output, nuis.rot_x ./ maxRot, nuis.rot_y ./ maxRot, nuis.rot_z ./ maxRot, ...
    nuis.rot_x_derivative1 ./ maxRotD1, nuis.rot_y_derivative1 ./ maxRotD1, nuis.rot_z_derivative1 ./ maxRotD1, ...
    nuis.trans_x ./ maxTrans, nuis.trans_y ./ maxTrans, nuis.trans_z ./ maxTrans, ...
    nuis.trans_x_derivative1 ./ maxTransD1, nuis.trans_y_derivative1 ./ maxTransD1, nuis.trans_z_derivative1 ./ maxTransD1];
% Replace NaNs with 0s
output(isnan(output)) = 0;

% Or export the entire table, but you run into issues w num cols:
% output = table2array(nuis);
end