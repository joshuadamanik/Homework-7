function X_star = minimize(X, eps, k, f)
    quasi = quasi_newton_class(length(X));
    for n = 1:1000
        p = quasi.bgfs(X, eps, f);
        [Xa, Xb] = unimodal_interval(X, eps, p, f);
        X_star = fibonacci_search(Xa, Xb, k, f);

        err = norm(X_star - X);
        if (err < eps * 0.01)
            break;
        end
        X = X_star;
    end
end