function ts = dynSeries(varargin)% dynSeries(a,b,c,d)

%@info:
%! @deftypefn {Function File} {@var{ts} =} dynSeries (@var{a},@var{b},@var{c},@var{d})
%! @anchor{dynSeries}
%! @sp 1
%! Constructor for the Dynare time series class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! T*1 vector or T*N matrix of data.
%! @item b
%! Initial date. For Quaterly, Monthly or Weekly data, b must be a string. For yearly data or if the frequence is not
%! defined b must be an integer.
%! @item c
%! N*q array of characters. Names of the N time series.
%! @item d
%! N*p array of characters. TeX names of the N time series.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Dynare time series object.
%! @end table
%! @sp 1
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item data
%! Array of doubles (nobs*vobs).
%! @item nobs
%! Scalar integer, the number of observations.
%! @item vobs
%! Scalar integer, the number of variables.
%! @item name
%! Array of chars (nvobs*n), names of the variables.
%! @item tex
%! Array of chars (nvobs*n), tex names of the variables.
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item time
%! Array of integers (nobs*2). The first column defines the years associated to each observation. The second column,
%! depending on the frequency, indicates the week, month or quarter numbers. For yearly data or unspecified frequency
%! the second column is filled by ones.
%! @item init
%! Row vector of integers (1*2) indicating the year and the week, month or quarter of the first observation. @var{init}
%! is the first row of @var{time}.
%! @item last
%! Row vector of integers (1*2) indicating the year and the week, month or quarter of the last observation. @var{init}
%! is the first row of @var{time}.
%! @end table
%! @sp 1
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{@@dynTime/dynTime}, @ref{@@dynTime/setTime}, @ref{@@dynTime/setFreq} 
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

ts = struct;

ts.data = [];
ts.nobs = 0;
ts.vobs = 0;
ts.name = {};
ts.tex  = {};
ts.freq = [];
ts.init = dynDate();
ts.time = dynDates();

ts = class(ts,'dynSeries');

switch nargin
  case 0
    %  Create an empty dynSeries object.
    return
  case 1
    if isa(varargin{1},'dynDate')
        if isempty(varargin{1})
            error(['dynSeries:: ' inputname(1) ' (identified as a dynDate object) must be non empty!'])
        else
            % Create an empty dynSeries object with an initial date.
            ts.init = varargin{1};
            ts.freq = varargin{1}.freq;
        end
        return
    elseif ischar(varargin{1})
        % Create a dynSeries object loading data in a file (*.csv, *.m, *.mat).
        if check_file_extension(varargin{1},'m')
            [freq,init,data,varlist] = load_m_file_data(varargin{1});
        elseif check_file_extension(varargin{1},'mat')
            [freq,init,data,varlist] = load_mat_file_data(varargin{1});
        elseif check_file_extension(varargin{1},'csv')
            [freq,init,data,varlist] = load_csv_file_data(varargin{1});
        else
            error(['dynSeries:: I''m not able to load data from ' inputname(1) '!'])
        end
        ts.init = init;
        ts.freq = freq;
        ts.data = data;
        ts.name = varlist;
        ts.vobs = length(varlist);
        ts.nobs = size(data,1); 
    end
  case {2,4}
    a = varargin{1};
    b = varargin{2};
    if nargin<4
        d = [];
    else
        d = varargin{4};
    end
    if nargin<3
        c = [];
    else
        c = varargin{3};
    end
    % Get data, number of observations and number of variables.
    ts.data = a;
    ts.nobs = size(a,1);
    ts.vobs = size(a,2);
    % Get the first date and set the frequency.
    if ~isempty(b)
        if ischar(b)% Weekly, Monthly or Quaterly data.
            quaterly = findstr('Q',b);
            monthly  = findstr('M',b);
            weekly   = findstr('W',b);
            yearly   = findstr('Y',b);
            if ~isempty(quaterly)
                ts.freq = 4;
            end
            if ~isempty(monthly)
                ts.freq = 12;
            end
            if ~isempty(weekly)
                ts.freq = 52;
            end
            if ~isempty(yearly)
                ts.freq = 1;
            end
            if isempty(quaterly) && isempty(monthly) && isempty(weekly)  && isempty(yearly)
                error('dynSeries:: Using a string as a second input argument, I can only handle weekly (W), monthly (M), quaterly (Q) or yearly (Y) data!');
            end
            ts.init = dynDate(b);
        elseif isa(b,'dynDate') && ~isempty(b)
            ts.freq = b.freq;
            ts.init = b;
        elseif isnumeric(b) && isreal(b) && isint(b)
            ts.freq = 1;
            ts.init = dynDate(b);
        else
            error('dynSeries::dynSeries: Wrong calling sequence!');
        end
    else% If b is empty.
        ts.freq = 1;
        ts.init = dynDate(1);
    end
    % Get the names of the variables.
    if ~isempty(c)
        if ts.vobs==length(c)
            for i=1:ts.vobs
                ts.name = vertcat(ts.name, c(i) );
            end
        else
            error('dynSeries::dynSeries: The number of declared names does not match the number of variables!')
        end
    else
        for i=1:ts.vobs
            ts.name = vertcat(ts.name, {'--NA--'});
        end
    end
    if ~isempty(d)
        if ts.vobs==length(d)
            for i=1:ts.vobs
                ts.tex = vertcat(ts.tex, d(i));
            end
        else
            error('dynSeries::dynSeries: The number of declared tex names does not match the number of variables!')
        end
    else
        for i=1:ts.vobs
            ts.tex = vertcat(ts.tex, {'--NA--'});
        end
    end
  otherwise
    error('dynSeries::dynSeries: Can''t instantiate the class, wrong calling sequence!')
end

ts.time = ts.init:(ts.init+ts.nobs);

%@test:1
%$ % Test if we can instantiate an empty dynSeries object.
%$ try
%$     ts = dynSeries();
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(4,1);
%$
%$ try
%$     aa = dynDate('1938M11');
%$     ts = dynSeries(aa);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,12);
%$     t(3) = dyn_assert(ts.init.freq,12);
%$     t(4) = dyn_assert(ts.init.time,[1938, 11]);
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dynSeries('dynseries_test_data.m');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1994, 3]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,100);
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dynSeries('dynseries_test_data.mat');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1994, 3]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,100);
%$ end
%$
%$ T = all(t);
%@eof:4

%@test:5
%$ t = zeros(7,1);
%$
%$ try
%$     ts = dynSeries('dynseries_test_data.csv');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1990, 1]);
%$     t(5) = dyn_assert(ts.vobs,4);
%$     t(6) = dyn_assert(ts.nobs,4);
%$     t(7) = dyn_assert(ts.name,{'azert';'yuiop';'qsdfg';'jklm'}); 
%$ end
%$
%$ T = all(t);
%@eof:5

%@test:6
%$ t = zeros(7,1);
%$
%$ try
%$     ts = dynSeries(transpose(1:5),[]);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,1);
%$     t(3) = dyn_assert(ts.init.freq,1);
%$     t(4) = dyn_assert(ts.init.time,[1, 1]);
%$     t(5) = dyn_assert(ts.vobs,1);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'--NA--'}); 
%$ end
%$
%$ T = all(t);
%@eof:6

%@test:7
%$ t = zeros(7,1);
%$
%$ try
%$     ts = dynSeries(transpose(1:5),'1950Q1');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,1);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'--NA--'}); 
%$ end
%$
%$ T = all(t);
%@eof:7


%@test:8
%$ t = zeros(7,1);
%$
%$ try
%$     ts = dynSeries([transpose(1:5), transpose(6:10)],'1950Q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Output'; 'Consumption'}); 
%$ end
%$
%$ T = all(t);
%@eof:8

